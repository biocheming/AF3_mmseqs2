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
        binary_path: str,
        database_path: str,
        n_cpu: int = 8,
        e_value: float = 0.0001,
        max_sequences: int = 10000,
        sensitivity: float = 7.5,
        gpu_devices: tuple[str, ...] | None = None,
    ):
        self.gpu_db = None  # <--- ADD THIS LINE HERE
        """Initializes the MMseqs2 runner.

        Args:
            binary_path: Path to the MMseqs2 binary.
            database_path: Path to the sequence database.
            n_cpu: Number of CPUs to use.
            e_value: E-value cutoff.
            max_sequences: Maximum number of sequences to return.
            sensitivity: Search sensitivity (from -s parameter).
            gpu_devices: List of GPU devices to use (e.g. ["0", "1"]).
        """
        self.binary_path = binary_path
        self.database_path = database_path
        self.n_cpu = n_cpu
        self.e_value = e_value
        self.max_sequences = max_sequences
        self.sensitivity = sensitivity
        self.use_gpu = gpu_devices is not None and len(gpu_devices) > 0
        self.gpu_devices = gpu_devices

        # 检查数据库文件是否存在
        if not os.path.exists(self.database_path):
            raise ValueError(f"Database file not found: {self.database_path}")

        # 获取数据库所在目录和基础名称
        db_dir = os.path.dirname(self.database_path)
        db_basename = os.path.basename(self.database_path)
        db_name = os.path.splitext(db_basename)[0]
        
        # 构建 MMseqs2 数据库目录结构
        self.mmseqs_dir = os.path.join(db_dir, f"{db_name}_mmseqs2")
        os.makedirs(self.mmseqs_dir, exist_ok=True)
        
        # 基础数据库文件
        self.base_db = os.path.join(self.mmseqs_dir, db_name)
        self.base_db_type = self.base_db + ".dbtype"
        self.base_db_source = self.base_db + ".source"

        # 如果基础数据库不存在，创建它
        if not os.path.exists(self.base_db_type):
            logging.info(f"Creating base MMseqs2 database at {self.base_db}")
            cmd = [
                self.binary_path,
                'createdb',
                self.database_path,
                self.base_db,
                '--dbtype', '1',  # 1表示氨基酸序列
                '--compressed', '0'
            ]
            try:
                run_with_logging(cmd)
                # 创建到原始文件的符号链接
                if not os.path.exists(self.base_db_source):
                    os.symlink(self.database_path, self.base_db_source)
            except subprocess.CalledProcessError as e:
                raise RuntimeError(f"Failed to create base database: {e}")

        # 如果启用了GPU，开始创建GPU索引
        if self.use_gpu:
            # GPU索引目录
            self.gpu_index_dir = os.path.join(self.mmseqs_dir, f"{db_name}_gpu_index")
            os.makedirs(self.gpu_index_dir, exist_ok=True)
            
            # GPU数据库文件
            self.gpu_db = os.path.join(self.gpu_index_dir, f"{db_name}_gpu")
            
            # 检查是否需要创建GPU索引
            required_exts = [".dbtype", ".index", ".lookup", "_h"]
            if not all(os.path.exists(self.gpu_db + ext) for ext in required_exts):
                logging.info(f"Creating GPU-optimized database in {self.gpu_index_dir}")
                
                # 1. 创建GPU优化的数据库
                cmd = [
                    self.binary_path,
                    'makepaddedseqdb',
                    self.base_db,
                    self.gpu_db
                ]
                try:
                    run_with_logging(cmd)
                except subprocess.CalledProcessError as e:
                    logging.error(f"Failed to create GPU-optimized database: {e}")
                    self.use_gpu = False
                    return

                # 2. 创建GPU优化的索引
                with tempfile.TemporaryDirectory() as tmp_dir:
                    cmd = [
                        self.binary_path,
                        'createindex',
                        self.gpu_db,
                        tmp_dir,
                        '--remove-tmp-files', '1',
                        '--threads', str(self.n_cpu),
                        '--comp-bias-corr', '0',  # 对GPU搜索禁用偏差校正
                        '--search-type', '1',     # 1表示氨基酸搜索
                        '--mask', '0'             # 对GPU搜索禁用掩码
                    ]
                    try:
                        run_with_logging(cmd)
                        logging.info("GPU index creation completed successfully")
                    except subprocess.CalledProcessError as e:
                        logging.error(f"Failed to create GPU index: {e}")
                        self.use_gpu = False

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
        gpu_db = self.gpu_db
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
                "--db-load-mode", "2",
                "--comp-bias-corr", "0",
                "--mask", "0",
                "--exact-kmer-matching", "1",
                "--use-gpu", "1"
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
                "--db-load-mode", "2"
            ]
            run_with_logging(cmd, env=env)
            
            return self._process_search_results(result_m8, target_sequence)
            
        finally:
            # 恢复原始环境变量
            if original_cuda_devices is not None:
                env["CUDA_VISIBLE_DEVICES"] = original_cuda_devices
            elif "CUDA_VISIBLE_DEVICES" in env:
                del env["CUDA_VISIBLE_DEVICES"]

    def _process_search_results(
        self,
        result_m8: str,
        target_sequence: str,
        query_db_path: str,
        tmp_dir: str,
    ) -> msa_tool.MsaToolResult:
        """Processes MMseqs2 search results from M8 to A3M format."""
        # Check if search produced results
        if not os.path.exists(result_m8) or os.path.getsize(result_m8) == 0:
            logging.error(f"No search results found in {result_m8}")
            return msa_tool.MsaToolResult(
                a3m="",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

        result_a3m_path = os.path.join(tmp_dir, "result.a3m")
        # In result2msa, the third argument is the target_db, which is self.base_db for CPU search
        # and gpu_db for GPU search. However, the prompt asks for self.base_db.
        # For now, sticking to self.base_db as per the prompt for the command structure.
        # The M8 file (result_m8) is the 4th argument to result2msa, not result_db.
        cmd = [
            self.binary_path,
            "result2msa",
            query_db_path, # Path to the query database
            self.base_db,   # Path to the target database (used for generating M8)
            result_m8,     # Path to the M8 file from convertalis
            result_a3m_path,
            "--db-load-mode", "2"
        ]
        try:
            run_with_logging(cmd)
            with open(result_a3m_path) as f:
                a3m_content = f.read()
            # Per prompt, only a3m content for success.
            # target_sequence and e_value are added for specific error cases.
            return msa_tool.MsaToolResult(a3m=a3m_content)
        except subprocess.CalledProcessError as e:
            logging.error(f"MMseqs2 result2msa failed: {e}")
            return msa_tool.MsaToolResult(a3m="")
        except IOError as e:
            logging.error(f"Failed to read A3M file {result_a3m_path}: {e}")
            return msa_tool.MsaToolResult(a3m="")

    def _cpu_search(
        self, 
        query_path: str, 
        result_m8: str, 
        tmp_dir: str,
        target_sequence: str,
    ) -> msa_tool.MsaToolResult:
        """Search using CPU-only MMseqs2."""
        query_db = os.path.join(tmp_dir, "query_db")
        result_db = os.path.join(tmp_dir, "result") # This is the raw search result, not M8
        
        # Create query DB
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        
        # Run CPU search
        cmd = [
            self.binary_path,
            "search",
            query_db,
            self.base_db,
            result_db, # Raw result database from search
            tmp_dir,
            "--threads", str(self.n_cpu),
            "-e", str(self.e_value),
            "--max-seqs", str(self.max_sequences),
            "-s", str(self.sensitivity),
            "--db-load-mode", "2"
        ]
        try:
            run_with_logging(cmd)
        except subprocess.CalledProcessError: # Errors in search itself
            logging.error("CPU search failed")
            # Return MsaToolResult with only a3m="" as per convention in existing code for search failure
            return msa_tool.MsaToolResult(a3m="")

        # Convert raw search results to m8 format
        # result_m8 is the path where the M8 file will be written
        cmd = [
            self.binary_path,
            "convertalis",
            query_db,
            self.base_db, # Target DB used in search
            result_db,   # Raw result DB from search
            result_m8,   # Output M8 file path
            "--format-mode", "0" # Output format m8
        ]
        try:
            run_with_logging(cmd)
        except subprocess.CalledProcessError: # Errors in M8 conversion
            logging.error("Failed to convert results to m8 format")
            return msa_tool.MsaToolResult(
                a3m="",
                target_sequence=target_sequence, # As per existing code for this error
                e_value=self.e_value
            )

        # Process the M8 file to A3M
        return self._process_search_results(
            result_m8=result_m8,
            target_sequence=target_sequence,
            query_db_path=query_db, # query_db is the path for the query database
            tmp_dir=tmp_dir
        )

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
            gpu_db = self.gpu_db
            if gpu_db:
                return self._gpu_search(query_path, result_m8, tmp_dir, target_sequence)
            else:
                logging.warning("GPU database indexing failed, falling back to CPU search")
            
            # 如果GPU搜索失败或未启用，使用CPU搜索
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
