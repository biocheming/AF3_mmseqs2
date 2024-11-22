"""MMseqs2 MSA tool."""

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Sequence

from alphafold3.data import parsers
from alphafold3.data.tools import msa_tool


def run_with_logging(cmd: Sequence[str], env: dict | None = None) -> None:
    """Runs command and logs stdout/stderr."""
    logging.info('Running command: %s', ' '.join(cmd))
    result = subprocess.run(
        cmd, capture_output=True, check=True, text=True,
        env=env
    )
    if result.stdout:
        logging.info("stdout:\n%s\n", result.stdout)
    if result.stderr:
        logging.info("stderr:\n%s\n", result.stderr)


class MMseqs2(msa_tool.MsaTool):
    """Python wrapper for MMseqs2."""

    def __init__(
        self,
        *,
        binary_path: str = "mmseqs",
        database_path: str,
        n_cpu: int = 8,
        e_value: float = 0.0001,
        max_sequences: int = 10000,
        sensitivity: float = 7.5,
        gpu_devices: list[str] | None = None,
    ):
        """Initialize MMseqs2 tool.
        
        Args:
            binary_path: Path to MMseqs2 binary
            database_path: Path to database FASTA file or MMseqs2 database
            n_cpu: Number of CPU threads to use
            e_value: E-value threshold for filtering hits
            max_sequences: Maximum number of sequences to return
            sensitivity: Search sensitivity (higher is more sensitive but slower)
            gpu_devices: List of GPU device indices to use (e.g. ["0", "1"])
        """
        self.binary_path = binary_path
        self.database_path = database_path
        self.n_cpu = n_cpu
        self.e_value = e_value
        self.max_sequences = max_sequences
        self.sensitivity = sensitivity
        self.gpu_devices = gpu_devices
        
        # 检查是否使用GPU
        self.use_gpu = bool(gpu_devices)
        if self.use_gpu:
            try:
                nvidia_smi = "nvidia-smi --query-gpu=gpu_bus_id --format=csv,noheader"
                result = subprocess.run(
                    nvidia_smi.split(),
                    capture_output=True,
                    check=True,
                    text=True
                )
                if result.stdout.strip():
                    logging.info(f"Found GPU devices: {result.stdout.strip()}")
                    logging.info(f"Will use GPU devices: {self.gpu_devices}")
                else:
                    logging.warning("No GPU devices found, falling back to CPU mode")
                    self.use_gpu = False
                    self.gpu_devices = None
            except (subprocess.SubprocessError, FileNotFoundError):
                logging.warning("Failed to detect GPU devices, falling back to CPU mode")
                self.use_gpu = False
                self.gpu_devices = None
        
        # 创建基础MMseqs2数据库
        if os.path.isfile(database_path):
            if database_path.endswith(('.fa', '.fasta')):
                self.database_path = self._create_base_db(database_path)
            else:
                # 检查是否是有效的MMseqs2数据库
                if not os.path.exists(database_path + ".index"):
                    logging.info(f"Creating MMseqs2 database from {database_path}")
                    self.database_path = self._create_base_db(database_path)
        else:
            # 尝试查找同名的FASTA文件
            fasta_path = database_path + ".fasta"
            if not os.path.exists(fasta_path):
                fasta_path = database_path + ".fa"
            
            if os.path.exists(fasta_path):
                logging.info(f"Found FASTA file at {fasta_path}")
                self.database_path = self._create_base_db(fasta_path)
            else:
                raise ValueError(f"Database not found at {database_path} or {database_path}.fasta")

    def _create_base_db(self, input_path: str) -> str:
        """创建基础的MMseqs2数据库（不包含索引）。

        Args:
            input_path: 输入文件路径

        Returns:
            MMseqs2数据库路径
        """
        # 获取输入文件的目录和文件名（不含扩展名）
        input_dir = os.path.dirname(input_path)
        input_basename = os.path.basename(input_path)
        db_name = os.path.splitext(input_basename)[0]
        
        # 创建数据库专属目录
        db_dir = os.path.join(input_dir, f"{db_name}_mmseqs2")
        os.makedirs(db_dir, exist_ok=True)
        
        # 设置数据库路径
        db_path = os.path.join(db_dir, db_name)
        
        # 如果数据库不存在则创建
        if not os.path.exists(db_path):
            logging.info(f"Creating base MMseqs2 database at {db_path}")
            try:
                cmd = [
                    self.binary_path,
                    "createdb",
                    input_path,
                    db_path,
                    '--dbtype', '1',
                    '--compressed', '0'
                ]
                run_with_logging(cmd)
            except subprocess.CalledProcessError as e:
                logging.error(f"Failed to create MMseqs2 database: {e}")
                raise
        return db_path

    def _ensure_database_indexed(self) -> str:
        """确保CPU版本的数据库索引存在。

        Returns:
            索引数据库的路径。
        """
        # 只在使用CPU模式时创建CPU索引
        if self.use_gpu:
            return self.database_path

        # 获取数据库所在目录
        db_dir = os.path.dirname(self.database_path)
        db_basename = os.path.basename(self.database_path)
        
        # 使用不带扩展名的数据库名作为索引前缀
        indexed_db = os.path.join(db_dir, "cpu_index", db_basename)
        os.makedirs(os.path.dirname(indexed_db), exist_ok=True)

        if not os.path.exists(indexed_db + ".index"):
            logging.info(f"Creating CPU database index in {os.path.dirname(indexed_db)}")
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                # 创建序列索引
                cmd = [
                    self.binary_path,
                    'createindex',
                    self.database_path,
                    tmp_dir,
                    '--remove-tmp-files', '1',
                    '--threads', str(self.n_cpu)
                ]
                run_with_logging(cmd)
        
        return indexed_db

    def _ensure_gpu_database_indexed(self) -> str | None:
        """确保GPU优化的数据库索引存在。

        Returns:
            GPU优化索引的路径, 如果创建失败则返回None。
        """
        if not self.use_gpu:
            return None

        # 获取数据库所在目录
        db_dir = os.path.dirname(self.database_path)
        db_basename = os.path.basename(self.database_path)
        
        # 使用不带扩展名的数据库名作为索引前缀
        gpu_index_dir = os.path.join(db_dir, "gpu_index")
        os.makedirs(gpu_index_dir, exist_ok=True)
        indexed_db = os.path.join(gpu_index_dir, db_basename)

        # 检查是否需要重新创建索引
        force_reindex = False
        if os.path.exists(indexed_db + ".idx"):
            try:
                verify_cmd = [
                    self.binary_path,
                    "touchdb",
                    indexed_db,
                    "--threads", "1"
                ]
                run_with_logging(verify_cmd)
            except subprocess.CalledProcessError:
                logging.warning("Existing GPU index appears to be corrupted, recreating...")
                force_reindex = True
                for ext in [".idx", ".idx.index", ".idx.dbtype"]:
                    if os.path.exists(indexed_db + ext):
                        os.remove(indexed_db + ext)
        
        if force_reindex or not all(os.path.exists(indexed_db + ext) 
                                  for ext in [".idx", ".idx.index", ".idx.dbtype"]):
            logging.info(f"Creating GPU-optimized database in {gpu_index_dir}")
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                # 1. 创建GPU优化的数据库
                cmd = [
                    self.binary_path,
                    'makepaddedseqdb',
                    self.database_path,
                    indexed_db
                ]
                run_with_logging(cmd)
                
                # 2. 创建GPU优化的索引
                cmd = [
                    self.binary_path,
                    'createindex',
                    indexed_db,
                    tmp_dir,
                    '--remove-tmp-files', '1',
                    '--threads', str(self.n_cpu),
                    '--comp-bias-corr', '0'
                ]
                
                try:
                    run_with_logging(cmd)
                    logging.info("GPU index creation completed successfully")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Failed to create GPU index: {e}")
                    return None
                
                # 3. 验证新创建的索引
                try:
                    verify_cmd = [
                        self.binary_path,
                        "touchdb",
                        indexed_db,
                        "--threads", "1"
                    ]
                    run_with_logging(verify_cmd)
                except subprocess.CalledProcessError:
                    logging.error("Failed to verify newly created GPU index")
                    return None
        
        return indexed_db

    def _get_source_db_path(self) -> str:
        """获取用于MSA转换的源数据库路径。"""
        db_basename = os.path.basename(self.database_path)
        db_name = os.path.splitext(db_basename)[0].replace('.fa', '').replace('.fasta', '')
        index_dir = os.path.join(os.path.dirname(self.database_path), f"{db_name}_mmseqs2_index")
        return os.path.join(index_dir, "source.fasta")

    def _gpu_search(
        self, 
        query_path: str, 
        result_m8: str, 
        tmp_dir: str,
        target_sequence: str,
    ) -> msa_tool.MsaToolResult:
        """Search using GPU-enabled MMseqs2."""
        # 获取GPU数据库路径
        gpu_db = self._ensure_gpu_database_indexed()
        if not gpu_db:
            logging.warning("GPU database not available, falling back to CPU search")
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
        
        query_db = os.path.join(tmp_dir, "query_db")
        result_db = os.path.join(tmp_dir, "result")
        
        # Create query DB
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        
        # Run GPU search with temporary environment variables
        env = os.environ.copy()
        original_cuda_devices = env.get("CUDA_VISIBLE_DEVICES")
        
        if self.gpu_devices:
            env["CUDA_VISIBLE_DEVICES"] = ",".join(self.gpu_devices)
            logging.info(f"Setting CUDA_VISIBLE_DEVICES={env['CUDA_VISIBLE_DEVICES']} for MMseqs2")
        
        try:
            cmd = [
                self.binary_path,
                "search",
                query_db,
                gpu_db,  # 使用GPU优化的数据库
                result_db,
                tmp_dir,
                "--threads", str(self.n_cpu),
                "--max-seqs", str(self.max_sequences),
                "-s", str(self.sensitivity),
                "-e", str(self.e_value),
                "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                "--format-mode", "2",
                "--db-load-mode", "2",
                "--comp-bias-corr", "0",
                "--mask", "0",
                "--orf-start-mode", "1",
                "--exact-kmer-matching", "1"
            ]
            run_with_logging(cmd, env=env)
            
            # Convert to m8 format
            cmd = [
                self.binary_path,
                "convertalis",
                query_db,
                gpu_db,
                result_db,
                result_m8,
                "--format-output", "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
                "--db-load-mode", "2",
                "--format-mode", "2"
            ]
            run_with_logging(cmd, env=env)
            
            return self._process_search_results(result_m8, target_sequence)
            
        finally:
            # 恢复原始环境变量
            if original_cuda_devices is not None:
                env["CUDA_VISIBLE_DEVICES"] = original_cuda_devices
            elif "CUDA_VISIBLE_DEVICES" in env:
                del env["CUDA_VISIBLE_DEVICES"]

    def _cpu_search(
        self, 
        query_path: str, 
        result_m8: str, 
        tmp_dir: str,
        target_sequence: str,
    ) -> msa_tool.MsaToolResult:
        """Search using CPU-only MMseqs2."""
        query_db = os.path.join(tmp_dir, "query_db")
        result_db = os.path.join(tmp_dir, "result")
        
        # Create query DB
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        
        # Run CPU search
        cmd = [
            self.binary_path,
            "search",
            query_db,
            self.database_path,
            result_db,
            tmp_dir,
            "--threads", str(self.n_cpu),
            "-e", str(self.e_value),
            "--max-seqs", str(self.max_sequences),
            "-s", str(self.sensitivity),
            "--db-load-mode", "2"
        ]
        try:
            run_with_logging(cmd)
        except subprocess.CalledProcessError:
            logging.error("CPU search failed")
            return msa_tool.MsaToolResult(a3m="")

        # Convert results to m8 format
        cmd = [
            self.binary_path,
            "convertalis",
            query_db,
            self.database_path,
            result_db,
            result_m8,
            "--format-mode", "0"
        ]
        try:
            run_with_logging(cmd)
        except subprocess.CalledProcessError:
            logging.error("Failed to convert results to m8 format")
            return msa_tool.MsaToolResult(
                a3m="",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

        # Check if search produced results
        if not os.path.exists(result_m8) or os.path.getsize(result_m8) == 0:
            logging.error("No search results found")
            return msa_tool.MsaToolResult(
                a3m="",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

        # Convert to a3m format
        result_a3m = os.path.join(tmp_dir, "result.a3m")
        cmd = [
            self.binary_path,
            "result2msa",
            query_db,
            self.database_path,
            result_db,
            result_a3m,
            "--db-load-mode", "2"
        ]
        try:
            run_with_logging(cmd)
            with open(result_a3m) as f:
                a3m = f.read()
            return msa_tool.MsaToolResult(a3m=a3m)
        except (subprocess.CalledProcessError, IOError) as e:
            logging.error(f"MSA conversion failed: {e}")
            return msa_tool.MsaToolResult(a3m="")

    def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
        """使用MMseqs2搜索序列数据库。

        Args:
            target_sequence: 目标序列。

        Returns:
            包含比对结果的MsaToolResult对象。
        """
        logging.info('Query sequence: %s', target_sequence)
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            query_path = os.path.join(tmp_dir, "query.fasta")
            result_m8 = os.path.join(tmp_dir, "result.m8")
            
            # 写入查询序列
            with open(query_path, "w") as f:
                f.write(f">query\n{target_sequence}\n")
            
            # 检测可用的GPU
            try:
                nvidia_smi = "nvidia-smi --query-gpu=gpu_bus_id --format=csv,noheader"
                result = subprocess.run(
                    nvidia_smi.split(),
                    capture_output=True,
                    check=True,
                    text=True
                )
                gpu_ids = result.stdout.strip().split('\n')
                if not gpu_ids:
                    logging.warning("No GPU devices found")
                    return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
                logging.info(f"Found {len(gpu_ids)} GPU devices")
            except (subprocess.SubprocessError, FileNotFoundError):
                logging.warning("Failed to get GPU information")
                return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

            # 如果使用GPU，确保有GPU索引
            gpu_db = self._ensure_gpu_database_indexed()
            if gpu_db:
                return self._gpu_search(query_path, result_m8, tmp_dir, target_sequence)
            else:
                logging.warning("GPU database indexing failed, falling back to CPU search")
            
            # 如果GPU搜索失败或未启用，使用CPU搜索
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
