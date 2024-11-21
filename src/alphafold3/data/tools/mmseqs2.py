"""MMseqs2 MSA tool."""

import logging
import os
import shutil
import subprocess
import tempfile
from typing import Sequence

from alphafold3.data import parsers
from alphafold3.data.tools import msa_tool


def run_with_logging(cmd: Sequence[str]) -> None:
    """Runs command and logs stdout/stderr."""
    logging.info('Running command: %s', ' '.join(cmd))
    result = subprocess.run(
        cmd, capture_output=True, check=True, text=True
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
        use_gpu: bool = False,
        e_value: float = 1e-4,
        max_sequences: int = 10_000,
        sensitivity: float = 7.5,
    ):
        """Initializes the wrapper.

        Args:
            binary_path: Path to the MMseqs2 binary.
            database_path: Path to the sequence database.
            n_cpu: Number of CPUs to use.
            use_gpu: Whether to use GPU acceleration.
            e_value: E-value threshold.
            max_sequences: Maximum number of sequences to return.
            sensitivity: Search sensitivity (from 1 to 7.5).
        """
        self.binary_path = binary_path
        self.database_path = database_path
        self.n_cpu = n_cpu
        self.use_gpu = use_gpu
        self.e_value = e_value
        self.max_sequences = max_sequences
        self.sensitivity = sensitivity

        # 检测可用的GPU
        if self.use_gpu:
            self.gpu_devices = self._get_gpu_devices()
            if not self.gpu_devices:
                logging.warning("No GPU devices found, falling back to CPU mode")
                self.use_gpu = False

    def _get_gpu_devices(self) -> list[int]:
        """获取可用的GPU设备ID列表。"""
        try:
            nvidia_smi = "nvidia-smi --query-gpu=gpu_bus_id --format=csv,noheader"
            result = subprocess.run(
                nvidia_smi.split(),
                capture_output=True,
                check=True,
                text=True
            )
            gpu_ids = list(range(len(result.stdout.strip().split('\n'))))
            logging.info(f"Found {len(gpu_ids)} GPU devices: {gpu_ids}")
            return gpu_ids
        except (subprocess.SubprocessError, FileNotFoundError):
            logging.warning("Failed to get GPU information")
            return []

    def _ensure_database_indexed(self) -> str:
        """确保数据库已建立索引。

        Returns:
            索引数据库的路径。
        """
        # 获取原始数据库文件名（不含路径）和扩展名
        db_basename = os.path.basename(self.database_path)
        db_name = os.path.splitext(db_basename)[0]
        
        # 创建该数据库专属的索引目录，移除任何.fa或.fasta后缀
        db_name = db_name.replace('.fa', '').replace('.fasta', '')
        index_dir = os.path.join(os.path.dirname(self.database_path), f"{db_name}_mmseqs2_index")
        os.makedirs(index_dir, exist_ok=True)
        
        indexed_db = os.path.join(index_dir, db_basename)
        
        if not os.path.exists(indexed_db) or not os.path.exists(indexed_db + ".index"):
            logging.info(f"Creating MMseqs2 database index in {index_dir}")
            
            # 创建到原始数据库的符号链接
            if os.path.exists(indexed_db):
                os.remove(indexed_db)
            os.symlink(self.database_path, indexed_db)
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                # 1. 创建标准MMseqs2数据库
                cmd = [
                    self.binary_path,
                    'createdb',
                    indexed_db,
                    indexed_db,
                    '--dbtype', '1',
                    '--compressed', '0'
                ]
                run_with_logging(cmd)
                
                # 2. 创建序列索引
                cmd = [
                    self.binary_path,
                    'createindex',
                    indexed_db,
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

        # 获取原始数据库文件名（不含路径）和扩展名
        db_basename = os.path.basename(self.database_path)
        db_name = os.path.splitext(db_basename)[0]
        
        # 使用该数据库专属的索引目录，移除任何.fa或.fasta后缀
        db_name = db_name.replace('.fa', '').replace('.fasta', '')
        gpu_index_dir = os.path.join(os.path.dirname(self.database_path), f"{db_name}_mmseqs2_index")
        os.makedirs(gpu_index_dir, exist_ok=True)
        
        gpu_index_base = os.path.join(gpu_index_dir, db_basename)
        
        # 检查是否需要重新创建索引
        force_reindex = False
        if os.path.exists(gpu_index_base + ".idx"):
            try:
                verify_cmd = [
                    self.binary_path,
                    "touchdb",
                    gpu_index_base,
                    "--threads", "1"
                ]
                run_with_logging(verify_cmd)
            except subprocess.CalledProcessError:
                logging.warning("Existing GPU index appears to be corrupted, recreating...")
                force_reindex = True
                for ext in [".idx", ".idx.index", ".idx.dbtype"]:
                    if os.path.exists(gpu_index_base + ext):
                        os.remove(gpu_index_base + ext)
        
        if force_reindex or not all(os.path.exists(gpu_index_base + ext) 
                                  for ext in [".idx", ".idx.index", ".idx.dbtype"]):
            logging.info(f"Creating GPU-optimized index in {gpu_index_dir}")
            
            with tempfile.TemporaryDirectory() as tmp_dir:
                # 根据数据库大小动态计算分块数
                db_size_gb = os.path.getsize(self._ensure_database_indexed()) / (1024**3)
                n_splits = max(20, int(db_size_gb / 2.5))  # 每2.5GB数据一个分块，最少20个分块
                
                index_cmd = [
                    self.binary_path,
                    "createindex",
                    self._ensure_database_indexed(),  # 使用已索引的数据库
                    gpu_index_base,
                    "--remove-tmp-files", "1",
                    "--threads", str(self.n_cpu),
                    "--comp-bias-corr", "0",
                    "--split", str(n_splits),
                    "--split-memory-limit", "2G"  # 每个split 2GB
                ]
                
                try:
                    run_with_logging(index_cmd)
                    logging.info("GPU index creation completed successfully")
                except subprocess.CalledProcessError as e:
                    logging.error(f"Failed to create GPU index: {e}")
                    return None
                
                # 验证新创建的索引
                try:
                    verify_cmd = [
                        self.binary_path,
                        "touchdb",
                        gpu_index_base,
                        "--threads", "1"
                    ]
                    run_with_logging(verify_cmd)
                except subprocess.CalledProcessError:
                    logging.error("Failed to verify newly created GPU index")
                    return None
        
        return gpu_index_base

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
            
            # 如果使用GPU，确保有GPU索引
            if self.use_gpu:
                gpu_db = self._ensure_gpu_database_indexed()
                if not gpu_db:
                    logging.warning("Failed to create/verify GPU index, falling back to CPU search")
                    return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
                
                success = False
                try:
                    gpu_cmd = [
                        self.binary_path,
                        'easy-search',
                        query_path,
                        gpu_db,
                        result_m8,
                        tmp_dir,
                        '--format-mode', '0',
                        '--threads', str(self.n_cpu),
                        '-e', str(self.e_value),
                        '--gpu', '1',
                        '--gpu-device', ','.join(map(str, self.gpu_devices)),
                        '--max-rejected', '5000',
                        '--max-seqs', str(self.max_sequences),
                        '-s', str(self.sensitivity),
                        '--remove-tmp-files', '1',
                        '--db-load-mode', '2'
                    ]
                    logging.info(f"Running GPU search with database: {gpu_db}")
                    run_with_logging(gpu_cmd)
                    success = True
                except subprocess.CalledProcessError as e:
                    logging.warning("GPU search failed: %s", e)
                
                if not success:
                    try:
                        gpu_cmd = [
                            self.binary_path,
                            'easy-search',
                            query_path,
                            gpu_db,
                            result_m8,
                            tmp_dir,
                            '--format-mode', '0',
                            '--threads', str(self.n_cpu),
                            '-e', str(self.e_value),
                            '--gpu', '1',
                            '--gpu-device', ','.join(map(str, self.gpu_devices)),
                            '--max-rejected', '3000',
                            '--max-seqs', '3000',
                            '-s', str(self.sensitivity),
                            '--remove-tmp-files', '1',
                            '--db-load-mode', '2'
                        ]
                        logging.info("Retrying GPU search with more aggressive memory settings")
                        run_with_logging(gpu_cmd)
                        success = True
                    except subprocess.CalledProcessError as e:
                        logging.warning("Second GPU search attempt failed: %s", e)
            
            if not self.use_gpu or not success:
                # 回退到CPU搜索
                return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
            
            # 将结果转换为a3m格式
            result_a3m = os.path.join(tmp_dir, "result.a3m")
            cmd = [
                self.binary_path,
                'result2msa',
                query_path,
                self._ensure_database_indexed(),
                result_m8,
                result_a3m,
                '--db-load-mode', '2'
            ]
            run_with_logging(cmd)
            
            with open(result_a3m) as f:
                a3m_content = f.read()
                return msa_tool.MsaToolResult(
                    target_sequence=target_sequence,
                    e_value=self.e_value,
                    a3m=a3m_content,
                )

    def _cpu_search(
        self, 
        query_path: str, 
        result_m8: str, 
        tmp_dir: str,
        target_sequence: str,
    ) -> msa_tool.MsaToolResult:
        """使用CPU模式进行搜索。"""
        logging.info("Running CPU-only search")
        cmd = [
            self.binary_path,
            'easy-search',
            query_path,
            self._ensure_database_indexed(),
            result_m8,
            tmp_dir,
            '--format-mode', '0',
            '--threads', str(self.n_cpu),
            '-e', str(self.e_value),
            '--max-seqs', str(self.max_sequences),
            '-s', str(self.sensitivity),
            '--split', '1',
            '--split-memory-limit', '30G',
        ]
        run_with_logging(cmd)
        
        # 转换为a3m格式
        result_a3m = os.path.join(tmp_dir, "result.a3m")
        cmd = [
            self.binary_path,
            'result2msa',
            query_path,
            self._ensure_database_indexed(),
            result_m8,
            result_a3m,
            '--db-load-mode', '2'
        ]
        run_with_logging(cmd)
        
        with open(result_a3m) as f:
            a3m_content = f.read()
            return msa_tool.MsaToolResult(
                target_sequence=target_sequence,
                e_value=self.e_value,
                a3m=a3m_content,
            )
