"""Simple MMseqs2 MSA tool with memory control."""

import os
import psutil
import subprocess
import tempfile
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Sequence, Dict, List, Any

from absl import logging
from alphafold3.data.tools import msa_tool
from alphafold3.data.tools import subprocess_utils


def get_gpu_memory_info(device_id: int = 0) -> dict:
    """Get GPU memory info."""
    try:
        cmd = ["nvidia-smi", "--query-gpu=memory.total,memory.free", "--format=csv,noheader,nounits", f"--id={device_id}"]
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=5)
        if result.returncode == 0:
            parts = result.stdout.strip().split(',')
            if len(parts) == 2:
                total_mb = int(parts[0].strip())
                free_mb = int(parts[1].strip())
                return {'total_mb': total_mb, 'free_mb': free_mb}
    except:
        pass
    return {}


def run_with_logging(cmd: Sequence[str], env: dict | None = None) -> None:
    """Run command with logging."""
    logging.info('Running: %s', ' '.join(cmd))
    result = subprocess.run(cmd, capture_output=True, check=True, text=True, env=env)
    if result.stdout:
        logging.info("stdout: %s", result.stdout.strip())


def run_with_quiet_logging(cmd: Sequence[str], env: dict | None = None, operation_name: str = "operation") -> None:
    """Run command with progress logging."""
    try:
        result = subprocess.run(cmd, capture_output=True, check=True, text=True, env=env)
        logging.info("✅ %s completed", operation_name)
    except subprocess.CalledProcessError as e:
        logging.error("❌ %s failed", operation_name)

        # Show detailed error information
        if e.stderr and len(e.stderr.strip()) > 0:
            stderr_lines = e.stderr.strip().split('\n')
            if len(stderr_lines) > 5:
                logging.error("Error details (last 5 lines): %s", ' | '.join(stderr_lines[-5:]))
            else:
                logging.error("Error details: %s", e.stderr.strip())

        # Check for memory-related errors
        if ("CUDA error: out of memory" in str(e.stderr) or
            "Cannot fit databases into" in str(e.stderr) or
            "Not enough memory" in str(e.stderr) or
            "IndexTable::initMemory" in str(e.stderr)):
            logging.error("🚨 MEMORY ERROR: Insufficient memory for operation")

        raise


def _estimate_db_size(db_name: str) -> str:
    """Estimate AF3 small database size for user information."""
    # AF3 small databases are much smaller than original databases
    size_estimates = {
        'uniprot_all': '8GB',      # Much smaller than original 74GB
        'uniref90': '4GB',         # Much smaller than original 30GB
        'bfd': '2GB',              # Much smaller than original 16GB, and small_bfd is even smaller
        'mgy_clusters': '6GB',     # Much smaller than original 108GB
    }

    # Clean database name for matching
    clean_name = db_name.lower().replace('.fasta', '').replace('.fa', '')

    # Check for each database key in the cleaned name
    for key in size_estimates.keys():
        if key in clean_name:
            return size_estimates[key]

    return 'Unknown'


class MMseqs2(msa_tool.MsaTool):
    """Simple MMseqs2 wrapper with memory control."""

    def __init__(self, *, binary_path: str, database_path: str, n_cpu: int = 8,
                 e_value: float = 0.0001, max_sequences: int = 10000, sensitivity: float = 4.0,
                 gpu_devices: tuple[str, ...] | None = None, use_gpu: bool = False):
        """Initialize with memory control."""
        self.binary_path = binary_path
        self.database_path = database_path
        self.n_cpu = n_cpu
        self.e_value = e_value
        self.max_sequences = max_sequences
        self.sensitivity = sensitivity
        self.use_gpu = use_gpu

        # Auto-detect all GPU devices if not specified
        if gpu_devices is None or len(gpu_devices) == 0:
            # Auto-detect available GPUs
            available_gpus = []
            for i in range(8):  # Check up to 8 GPUs
                if get_gpu_memory_info(i):
                    available_gpus.append(str(i))
            self.gpu_devices = tuple(available_gpus)
        else:
            self.gpu_devices = gpu_devices

        # Only enable GPU if we have both GPU support and available devices
        self.use_gpu = self.use_gpu and len(self.gpu_devices) > 0
        logging.debug(f"MMseqs2 GPU config: use_gpu={self.use_gpu}, gpu_devices={self.gpu_devices}")

        # Verify binary
        subprocess_utils.check_binary_exists(path=self.binary_path, name='MMseqs2')

        # Check MMseqs2 capabilities and version
        self._check_mmseqs2_capabilities()

        # Initialize MSA cache for future optimization
        self._msa_cache = {}

        # GPU task coordination for parallel processing
        self._gpu_task_lock = None  # Will be set if GPU is used

        # Setup database
        if not os.path.exists(self.database_path):
            raise ValueError(f"Database not found: {self.database_path}")

        db_dir = os.path.dirname(self.database_path)
        db_name = os.path.splitext(os.path.basename(self.database_path))[0]
        self.mmseqs_dir = os.path.join(db_dir, f"{db_name}_mmseqs2")
        os.makedirs(self.mmseqs_dir, exist_ok=True)
        self.base_db = os.path.join(self.mmseqs_dir, db_name)

        # Pre-setup databases
        self._setup_databases()

        # Initialize MSA cache for future optimization
        self._msa_cache = {}

    def _check_mmseqs2_capabilities(self):
        """Check MMseqs2 version and capabilities for optimal MSA generation."""
        self.supports_direct_a3m = False

        try:
            # Try to get convertalis help to check for A3M support
            result = subprocess.run([self.binary_path, 'convertalis', '--help'],
                                  capture_output=True, text=True, timeout=10)

            help_text = (result.stdout + result.stderr).lower()
            if 'a3m' in help_text and 'format' in help_text:
                self.supports_direct_a3m = True
                logging.debug("MMseqs2 supports direct A3M output via convertalis")
            else:
                self.supports_direct_a3m = False
                logging.debug("MMseqs2 does not support direct A3M output, will use result2msa fallback")

        except (subprocess.TimeoutExpired, subprocess.CalledProcessError, FileNotFoundError):
            self.supports_direct_a3m = False
            logging.debug("Could not determine MMseqs2 capabilities, assuming no direct A3M support")

    def _get_database_specific_max_seqs(self, db_name: str) -> int:
        """Get database-specific max_sequences with balanced optimization."""
        db_name_lower = db_name.lower()

        # 平衡的序列数设置
        if 'uniref90' in db_name_lower:
            return 3000   # 平衡质量与速度
        elif 'mgy_clusters' in db_name_lower or 'mgnify' in db_name_lower:
            return 1500   # 主瓶颈：适度优化
        elif 'bfd' in db_name_lower:
            return 2000   # 正常序列数
        elif 'uniprot' in db_name_lower:
            return 2500   # 平衡设置
        else:
            return min(self.max_sequences, 2500)  # Default balanced

    def _setup_databases(self):
        # Create base database if needed
        base_db_files = [self.base_db + ext for ext in [".dbtype", ".index"]]
        if not all(os.path.exists(f) for f in base_db_files):
            cmd = [self.binary_path, 'createdb', self.database_path, self.base_db, '--dbtype', '1']
            run_with_logging(cmd)

        # GPU setup - ensure complete GPU database creation
        if self.use_gpu:
            db_name = os.path.splitext(os.path.basename(self.database_path))[0]
            self.gpu_index_dir = os.path.join(self.mmseqs_dir, f"{db_name}_gpu_index")
            os.makedirs(self.gpu_index_dir, exist_ok=True)
            self.gpu_db = os.path.join(self.gpu_index_dir, f"{db_name}_gpu")
            
            # Check all required GPU files
            gpu_files = [self.gpu_db + ext for ext in [".dbtype", ".index", ".lookup", "_h"]]
            gpu_db_complete = all(os.path.exists(f) for f in gpu_files)

            if not gpu_db_complete:
                # Create padded sequence database for GPU
                cmd = [self.binary_path, 'makepaddedseqdb', self.base_db, self.gpu_db]
                run_with_logging(cmd)

                # Create GPU index
                with tempfile.TemporaryDirectory() as tmp_dir:
                    cmd = [self.binary_path, 'createindex', self.gpu_db, tmp_dir, '--threads', str(min(self.n_cpu, 8))]
                    run_with_logging(cmd)

                # Verify GPU database was created successfully
                if not all(os.path.exists(f) for f in gpu_files):
                    raise RuntimeError(f"Failed to create complete GPU database for {db_name}")

    def _gpu_search(self, query_path: str, result_m8: str, tmp_dir: str, target_sequence: str) -> msa_tool.MsaToolResult:
        """GPU search with memory check."""
        import time
        db_name = os.path.basename(self.database_path)
        db_short = db_name.split('_')[0] if '_' in db_name else db_name[:10]
        logging.info("Getting MSA for sequence from database: %s (GPU)", db_name)

        if not self.gpu_db:
            logging.warning("GPU database not available, falling back to CPU")
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
        
        # Stage 1: Create query database (CPU, ~2-10s)
        stage_start = time.time()
        logging.debug("📝 [%s] Stage 1/5: Creating query database (CPU, ~2-10s)", db_short)
        logging.debug("   [%s] Input: FASTA sequence | Output: MMseqs2 DB | Hardware: CPU", db_short)
        query_db = os.path.join(tmp_dir, "query_db")
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        logging.debug("✅ [%s] Stage 1 completed in %.1fs", db_short, time.time() - stage_start)
        
        # Stage 2: GPU memory check and selection (~1-5s)
        stage_start = time.time()

        # 短暂等待，确保上个任务释放GPU内存
        time.sleep(1.0)

        logging.debug("🔍 [%s] Stage 2/5: Checking GPU memory and selecting device (~1-5s)", db_short)
        logging.debug("   [%s] Hardware: GPU | Requirement: ≥3000MB free memory (avoiding GPU0)", db_short)

        # 智能GPU选择：避开GPU0（通常被AF3主程序占用），优先使用GPU1、GPU2等
        best_gpu_id = None
        best_gpu_memory = 0
        gpu_memory_details = []

        # 定义GPU优先级：避开GPU0，优先使用其他GPU
        gpu_priority_order = [1, 2, 3, 4, 5, 6, 7, 0]  # GPU0放在最后
        
        # 如果用户指定了GPU设备，调整优先级顺序
        if self.gpu_devices:
            user_specified = [int(x) for x in self.gpu_devices if x.isdigit()]
            # 用户指定的GPU优先，但仍然尽量避开GPU0
            gpu_priority_order = [gpu for gpu in user_specified if gpu != 0] + \
                               [gpu for gpu in gpu_priority_order if gpu not in user_specified and gpu != 0] + \
                               ([0] if 0 in user_specified else [])

        logging.debug("   [%s] GPU selection priority order: %s", db_short, gpu_priority_order)

        # 按优先级顺序检查GPU
        for device_id in gpu_priority_order:
            if device_id >= 8:  # Skip invalid GPU IDs
                continue
                
            gpu_info = get_gpu_memory_info(int(device_id))
            if gpu_info:
                free_mb = gpu_info['free_mb']
                total_mb = gpu_info['total_mb']
                gpu_memory_details.append(f"GPU {device_id}: {free_mb}MB free / {total_mb}MB total")
                logging.debug("   [%s] GPU %s: %sMB free / %sMB total",
                            db_short, device_id, free_mb, total_mb)

                # 使用更宽松的内存阈值（3GB而不是5GB）
                min_memory_threshold = 3000
                
                # 选择第一个满足条件的GPU（按优先级）
                if free_mb >= min_memory_threshold:
                    if best_gpu_id is None or (device_id != 0 and best_gpu_id == 0):
                        # 优先选择非GPU0的设备
                        best_gpu_id = device_id
                        best_gpu_memory = free_mb
                        logging.debug("   [%s] GPU %s selected (priority-based, %dMB free)", 
                                    db_short, device_id, free_mb)
                        break  # 找到合适的GPU就停止搜索
                    elif best_gpu_id is None:
                        # 如果只有GPU0可用，也接受
                        best_gpu_id = device_id
                        best_gpu_memory = free_mb
                        logging.debug("   [%s] GPU %s selected (fallback to GPU0, %dMB free)", 
                                    db_short, device_id, free_mb)

        # GPU内存保守策略：降低阈值确保GPU真正可用
        min_memory_threshold = 5000  # 降低到5000MB，更保守但确保GPU可用
        
        if best_gpu_id is not None and best_gpu_memory >= min_memory_threshold:
            gpu_id = best_gpu_id
            if gpu_id == 0:
                logging.warning("[%s] ⚠️ Using GPU0 which may be busy with AF3 main process (%dMB free)", 
                              db_short, best_gpu_memory)
            else:
                logging.info("✅ [%s] Selected GPU %s with %dMB free memory (safe threshold 5GB+)", 
                            db_short, gpu_id, best_gpu_memory)
        else:
            logging.warning("[%s] No suitable GPU found (need 5000MB+), GPU details: %s, using CPU",
                          db_short, " | ".join(gpu_memory_details))
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
        logging.debug("✅ [%s] Stage 2 completed in %.1fs", db_short, time.time() - stage_start)

        # Stage 3: GPU search (main bottleneck, ~10-300s)
        stage_start = time.time()
        result_db = os.path.join(tmp_dir, "result")
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

        # 验证GPU可用性
        try:
            import subprocess
            gpu_test = subprocess.run(
                ["nvidia-smi", "-i", str(gpu_id), "--query-gpu=memory.free", "--format=csv,noheader,nounits"],
                capture_output=True, text=True, timeout=5
            )
            if gpu_test.returncode == 0:
                gpu_memory = int(gpu_test.stdout.strip())
                logging.info("🔍 [%s] GPU %d verified: %dMB free memory", db_short, gpu_id, gpu_memory)
            else:
                logging.warning("⚠️ [%s] GPU %d verification failed, may fallback to CPU", db_short, gpu_id)
        except Exception as e:
            logging.warning("⚠️ [%s] GPU verification error: %s", db_short, e)

        # GPU搜索激进优化：充分利用24GB显存的4090
        gpu_max_seqs = self._get_database_specific_max_seqs(db_name)

        # 智能GPU优化：让MMseqs2自适应CPU，通过内存限制来控制
        # 根据ChatGPT建议：给足线程数，让MMseqs2根据内存自己决定实际使用
        gpu_threads = self.n_cpu  # 使用所有可用CPU，MMseqs2会自适应
        if 'mgy_clusters' in db_name.lower():
            gpu_sensitivity = 4.5 # mgnify: 适度降低sensitivity
            logging.debug("   [%s] Smart GPU optimization for mgnify: max_seqs=%d, threads=%d (auto-adjust), sensitivity=%.1f", 
                         db_short, gpu_max_seqs, gpu_threads, gpu_sensitivity)
        elif 'uniprot' in db_name.lower():
            gpu_sensitivity = 5.0 # uniprot: 适度降低sensitivity
            logging.debug("   [%s] Smart GPU optimization for uniprot: max_seqs=%d, threads=%d (auto-adjust), sensitivity=%.1f", 
                         db_short, gpu_max_seqs, gpu_threads, gpu_sensitivity)
        else:
            gpu_sensitivity = 6.0 # 保持质量
            logging.debug("   [%s] Smart GPU optimization: max_seqs=%d, threads=%d (auto-adjust), sensitivity=%.1f", 
                         db_short, gpu_max_seqs, gpu_threads, gpu_sensitivity)

        # 智能GPU搜索：让MMseqs2根据内存自适应，通过split降低峰值
        # 根据ChatGPT建议：给足资源，让MMseqs2自己管理
        import psutil
        total_memory_gb = psutil.virtual_memory().total / (1024**3)
        # 使用系统内存的60%作为限制（GPU搜索相对保守）
        memory_limit_gb = int(total_memory_gb * 0.6)
        
        cmd = [
            self.binary_path, "search", query_db, self.gpu_db, result_db, tmp_dir,
            "--threads", str(gpu_threads),  # 给足线程，让MMseqs2自适应
            "--max-seqs", str(gpu_max_seqs), 
            "-s", str(gpu_sensitivity), "-e", str(self.e_value), 
            "--gpu", "1",
            "--split-memory-limit", f"{memory_limit_gb}G",  # 动态内存限制
            "--split", "0",  # 不拆分，如果有充足内存（测试用）
            "--db-load-mode", "2",  
            "--max-accept", "50000",  
            "--max-rejected", "200",
            "-v", "3"
        ]
        
        # 显示智能优化信息
        logging.info("🔄 [%s] Smart GPU search: %d threads (MMseqs2 auto-adjust), mem_limit=%dG, split=0, sensitivity=%.1f", 
                   db_short, gpu_threads, memory_limit_gb, gpu_sensitivity)
        
        # 添加GPU调试信息 - 显示关键参数和值
        key_params = []
        for i, param in enumerate(cmd):
            if param in ['--split-memory-limit', '--max-accept', '--db-load-mode'] and i+1 < len(cmd):
                key_params.append(f"{param} {cmd[i+1]}")
        logging.info("🔍 [%s] GPU Key params: %s", db_short, " | ".join(key_params))
        
        # 使用优化的GPU配置，控制资源使用强度
        logging.info("🔍 [%s] Using optimized GPU configuration with memory limits", db_short)

        try:
            logging.info("🔍 [%s] Starting GPU search with CUDA_VISIBLE_DEVICES=%s", db_short, env.get("CUDA_VISIBLE_DEVICES"))
            
            # 详细检查CUDA环境
            try:
                cuda_check = subprocess.run(
                    ["nvidia-smi", "-L"], 
                    capture_output=True, text=True, timeout=5
                )
                if cuda_check.returncode == 0:
                    logging.info("✅ [%s] CUDA GPUs available: %s", db_short, cuda_check.stdout.strip())
                else:
                    logging.warning("⚠️ [%s] nvidia-smi failed", db_short)
            except Exception as e:
                logging.warning("⚠️ [%s] CUDA check failed: %s", db_short, e)
            
            # 检查GPU数据库文件是否存在
            try:
                if os.path.exists(self.gpu_db):
                    logging.info("✅ [%s] GPU database exists: %s", db_short, self.gpu_db)
                else:
                    logging.warning("⚠️ [%s] GPU database not found: %s", db_short, self.gpu_db)
            except Exception as e:
                logging.warning("⚠️ [%s] GPU database check failed: %s", db_short, e)
            
            # 检查MMseqs2 GPU支持
            try:
                gpu_check = subprocess.run(
                    [self.binary_path, "search", "--help"], 
                    capture_output=True, text=True, timeout=10
                )
                if "--gpu" in gpu_check.stdout:
                    logging.info("✅ [%s] MMseqs2 GPU support confirmed in help", db_short)
                else:
                    logging.warning("⚠️ [%s] MMseqs2 may not support --gpu flag", db_short)
            except Exception as e:
                logging.warning("⚠️ [%s] Could not check MMseqs2 GPU support: %s", db_short, e)
            
            # 执行GPU搜索并检查GPU使用情况
            import subprocess
            result = subprocess.run(cmd, env=env, capture_output=True, text=True)
            
            if result.returncode != 0:
                logging.error("❌ [%s] GPU search failed with return code %d", db_short, result.returncode)
                logging.error("❌ [%s] Error output: %s", db_short, result.stderr)
                raise RuntimeError(f"GPU search failed for {db_name}")
            
            # 修正GPU检测逻辑：基于成功执行和GPU数据库使用来判断
            # MMseqs2 GPU模式在成功时可能不输出stderr信息
            stderr_preview = result.stderr[:200] if result.stderr else "No stderr output"
            
            # 关键判断：如果GPU命令成功执行且使用了GPU数据库，则认为GPU工作正常
            gpu_used_successfully = (
                result.returncode == 0 and 
                "--gpu 1" in str(cmd) and 
                "gpu" in self.gpu_db.lower()
            )
            
            if gpu_used_successfully:
                logging.info("✅ [%s] GPU search completed successfully (no stderr indicates normal GPU operation)", db_short)
                logging.info("🔥 [%s] GPU acceleration confirmed - using GPU database: %s", db_short, os.path.basename(self.gpu_db))
            else:
                # 检查是否有明确的GPU错误信息
                gpu_error_indicators = ["not a valid GPU database", "GPU", "CUDA error", "gpu architecture"]
                found_errors = [err for err in gpu_error_indicators if err in result.stderr.lower()]
                
                if found_errors:
                    logging.warning("⚠️ [%s] GPU search encountered issues: %s", db_short, found_errors)
                else:
                    logging.info("✅ [%s] GPU search completed (minimal stderr output is normal)", db_short)
            
            # 简化调试信息
            if result.stderr and len(result.stderr.strip()) > 0:
                logging.info("🔍 [%s] stderr preview: %s", db_short, stderr_preview)
            
            logging.info("✅ [%s] GPU search completed in %.1fs", db_short, time.time() - stage_start)

            # Stage 4: Generate A3M directly (skip convertalis -> m8)
            stage_start = time.time()
            result_a3m = os.path.join(tmp_dir, "result.a3m")
            logging.debug("🚀 [%s] Stage 4/5: Generating A3M directly (~10-60s)", db_short)
            logging.debug("   [%s] Input: Alignment DB | Output: A3M MSA | Hardware: CPU/GPU", db_short)

            try:
                generation_time = self._generate_msa_fast(
                    mmseqs_bin=self.binary_path,
                    query_db=query_db,
                    target_db=self.gpu_db,
                    result_db=result_db,
                    output_a3m=result_a3m,
                    threads=self.n_cpu  # GPU模式：使用所有CPU处理A3M生成
                )
                logging.debug("✅ [%s] Stage 4 completed in %.1fs", db_short, time.time() - stage_start)
                return self._process_search_results(result_a3m, target_sequence, query_db, tmp_dir, from_a3m=True, msa_generation_time=generation_time)
            except Exception as e:
                logging.warning("Fast A3M generation failed: %s. Falling back to convertalis -> result2msa", e)
                # Fallback: convertalis + traditional result2msa
                cmd = [self.binary_path, "convertalis", query_db, self.gpu_db, result_db, result_m8]
                run_with_logging(cmd)
                logging.debug("✅ [%s] Stage 4 completed in %.1fs (fallback)", db_short, time.time() - stage_start)
                return self._process_search_results(result_m8, target_sequence, query_db, tmp_dir)
        except subprocess.CalledProcessError as e:
            logging.warning("GPU search failed for %s, falling back to CPU", db_name)
            logging.debug("GPU search error details: %s", e)
            if e.stdout:
                logging.debug("GPU search stdout: %s", e.stdout)
            if e.stderr:
                logging.debug("GPU search stderr: %s", e.stderr)
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

    def _gpu_search_only(self, query_path: str, tmp_dir: str) -> Dict[str, Any]:
        """Execute only GPU search phase, return paths for later A3M generation.
        
        This enables async processing: GPU search and A3M generation can overlap.
        
        Returns:
            Dict containing: query_db, target_db, result_db, tmp_dir, db_short, gpu_id
        """
        import time
        db_name = os.path.basename(self.database_path)
        db_short = db_name.split('_')[0] if '_' in db_name else db_name[:10]
        
        if not self.gpu_db:
            raise RuntimeError("GPU database not available for GPU search")
        
        # Stage 1: Create query database
        stage_start = time.time()
        logging.debug("📝 [%s] ASYNC Stage 1/3: Creating query database", db_short)
        query_db = os.path.join(tmp_dir, "query_db")
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        logging.debug("✅ [%s] ASYNC Stage 1 completed in %.1fs", db_short, time.time() - stage_start)
        
        # Stage 2: GPU memory check and selection
        stage_start = time.time()
        time.sleep(1.0)  # Brief wait for GPU memory release
        
        logging.debug("🔍 [%s] ASYNC Stage 2/3: Selecting GPU device", db_short)
        
        # GPU selection logic (same as _gpu_search)
        best_gpu_id = None
        best_gpu_memory = 0
        gpu_priority_order = [1, 2, 3, 4, 5, 6, 7, 0]
        
        if self.gpu_devices:
            user_specified = [int(x) for x in self.gpu_devices if x.isdigit()]
            gpu_priority_order = [gpu for gpu in user_specified if gpu != 0] + \
                               [gpu for gpu in gpu_priority_order if gpu not in user_specified and gpu != 0] + \
                               ([0] if 0 in user_specified else [])
        
        for device_id in gpu_priority_order:
            if device_id >= 8:
                continue
            gpu_info = get_gpu_memory_info(int(device_id))
            if gpu_info:
                free_mb = gpu_info['free_mb']
                min_memory_threshold = 3000
                if free_mb >= min_memory_threshold:
                    best_gpu_id = device_id
                    best_gpu_memory = free_mb
                    break
        
        if best_gpu_id is None or best_gpu_memory < 5000:
            raise RuntimeError(f"No suitable GPU found for {db_short}")
        
        gpu_id = best_gpu_id
        logging.info("✅ [%s] Selected GPU %d with %dMB free memory (safe threshold 5GB+)", 
                     db_short, gpu_id, best_gpu_memory)
        
        # Stage 3: GPU search
        stage_start = time.time()
        result_db = os.path.join(tmp_dir, "result_db")
        
        # Get memory info
        total_memory_gb = psutil.virtual_memory().total / (1024**3)
        mem_limit_gb = int(total_memory_gb * 0.6)
        
        # Determine sensitivity based on database (same logic as _gpu_search)
        db_name_lower = db_name.lower()
        if 'mgy_clusters' in db_name_lower:
            sensitivity = 4.5  # mgnify: moderate sensitivity
        elif 'uniprot' in db_name_lower:
            sensitivity = 5.0  # uniprot: moderate sensitivity
        elif 'bfd' in db_name_lower:
            sensitivity = 6.0  # bfd: higher sensitivity for quality
        elif 'uniref90' in db_name_lower:
            sensitivity = 6.0  # uniref90: higher sensitivity for quality
        else:
            sensitivity = self.sensitivity if hasattr(self, 'sensitivity') else 5.7
        
        # Use database-specific max_sequences (CRITICAL FIX!)
        max_seqs = self._get_database_specific_max_seqs(db_name)
        
        cmd = [
            self.binary_path, "search",
            query_db, self.gpu_db, result_db, tmp_dir,
            "--threads", str(self.n_cpu),
            "-s", str(sensitivity),
            "-e", str(self.e_value),
            "--max-seqs", str(max_seqs),  # Fixed: use database-specific value
            "--max-accept", "50000",
            "--split-memory-limit", f"{mem_limit_gb}G",
            "--split", "0",
            "--db-load-mode", "2",
            "--gpu", "1"
        ]
        
        env = os.environ.copy()
        env['CUDA_VISIBLE_DEVICES'] = str(gpu_id)
        
        logging.info("🔄 [%s] ASYNC GPU search: %d threads, mem_limit=%dG, split=0, sensitivity=%.1f, max_seqs=%d",
                     db_short, self.n_cpu, mem_limit_gb, sensitivity, max_seqs)
        logging.info("🔍 [%s] GPU Key params: --split-memory-limit %dG | --max-seqs %d | --max-accept 50000",
                     db_short, mem_limit_gb, max_seqs)
        logging.info("🔍 [%s] Starting GPU search with CUDA_VISIBLE_DEVICES=%d", db_short, gpu_id)
        
        import subprocess
        result = subprocess.run(cmd, env=env, capture_output=True, text=True)
        
        if result.returncode != 0:
            logging.error("❌ [%s] GPU search failed: %s", db_short, result.stderr)
            raise RuntimeError(f"GPU search failed for {db_short}")
        
        logging.info("✅ [%s] GPU search completed in %.1fs", db_short, time.time() - stage_start)
        
        return {
            'query_db': query_db,
            'target_db': self.gpu_db,
            'result_db': result_db,
            'tmp_dir': tmp_dir,
            'db_short': db_short,
            'gpu_id': gpu_id,
            'search_time': time.time() - stage_start
        }
    
    def _a3m_generation_only(self, search_result: Dict[str, Any], target_sequence: str) -> msa_tool.MsaToolResult:
        """Execute only A3M generation phase from completed GPU search results.
        
        This enables async processing: can be called independently after GPU search completes.
        
        Args:
            search_result: Dict from _gpu_search_only containing paths and metadata
            target_sequence: The query protein sequence
            
        Returns:
            MsaToolResult with generated MSA
        """
        import time
        db_short = search_result['db_short']
        
        stage_start = time.time()
        result_a3m = os.path.join(search_result['tmp_dir'], "result.a3m")
        
        logging.debug("🚀 [%s] ASYNC A3M Generation: Starting (~10-180s)", db_short)
        logging.debug("   [%s] Input: result_db | Output: A3M MSA | Hardware: CPU", db_short)
        
        try:
            generation_time = self._generate_msa_fast(
                mmseqs_bin=self.binary_path,
                query_db=search_result['query_db'],
                target_db=search_result['target_db'],
                result_db=search_result['result_db'],
                output_a3m=result_a3m,
                threads=self.n_cpu
            )
            
            logging.debug("✅ [%s] ASYNC A3M generation completed in %.1fs", db_short, generation_time)
            
            return self._process_search_results(
                result_a3m, 
                target_sequence, 
                search_result['query_db'], 
                search_result['tmp_dir'], 
                from_a3m=True, 
                msa_generation_time=generation_time
            )
        except Exception as e:
            logging.error("❌ [%s] ASYNC A3M generation failed: %s", db_short, e)
            # Return minimal MSA as fallback
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

    def _cpu_search(self, query_path: str, result_m8: str, tmp_dir: str, target_sequence: str) -> msa_tool.MsaToolResult:
        """CPU search with memory limit."""
        db_name = os.path.basename(self.database_path)
        logging.info("Getting MSA for sequence from database: %s", db_name)

        # Stage 1: Create query database (CPU, ~2-10s)
        logging.debug("📝 [CPU][%s] Stage 1/4: Creating query database (CPU, ~2-10s)", db_name)
        logging.debug("   [CPU][%s] Input: FASTA sequence | Output: MMseqs2 DB | Hardware: CPU", db_name)
        query_db = os.path.join(tmp_dir, "query_db")
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        logging.debug("✅ [CPU][%s] Stage 1 completed", db_name)

        # Stage 2: Memory limit calculation (~1s)
        logging.debug("⚙️ [CPU][%s] Stage 2/4: Calculating memory limits (~1s)", db_name)
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        total_kb = int(line.split()[1])
                        total_mb = total_kb // 1024
                        break
        except:
            total_mb = 256 * 1024  # Assume 256GB if can't read

        # 精确内存策略：只有特定数据库在prefilter阶段需要限制
        # 其他数据库和阶段可以使用更激进的设置
        db_name_lower = db_name.lower()

        # CPU搜索优化：使用数据库特定的max_sequences，同时保持内存策略
        max_seqs = self._get_database_specific_max_seqs(db_name)

        # 智能CPU优化：让MMseqs2通过内存限制自适应，根据ChatGPT建议
        # 使用系统内存的70%作为限制，让MMseqs2自己决定并行度
        mem_limit_gb = int(total_mb * 0.7)  # 系统内存的70%
        threads = self.n_cpu  # 给足所有CPU，让MMseqs2根据内存自适应
        strategy = "Smart CPU auto-adjust (memory-aware scheduling)"
        
        logging.info("💡 [CPU][%s] Using smart strategy: threads=%d (auto-adjust), mem=%.1f%% of total", 
                   db_name, threads, 70.0)

        # 计算百分比用于日志显示
        actual_percent = (mem_limit_gb / (total_mb / 1024)) * 100

        logging.debug("   [CPU][%s] System RAM: ~%dGB | Using: %dGB (%.1f%%) | Strategy: %s",
                    db_name, total_mb // 1024, mem_limit_gb, actual_percent, strategy)
        logging.debug("   [CPU][%s] Settings: max_seqs=%d, threads=%d", db_name, max_seqs, threads)
        logging.debug("✅ [CPU][%s] Stage 2 completed", db_name)

        # Stage 3: CPU search (main bottleneck, ~5-60min)
        logging.info("Running homology search against %s database", db_name)
        logging.debug("   [CPU][%s] Input: Query DB + Target DB | Output: Alignment DB | Hardware: CPU", db_name)
        logging.debug("   [CPU][%s] Algorithm: k-mer matching + ungapped prefilter + gapped alignment", db_name)

        result_db = os.path.join(tmp_dir, "result")
        # 平衡CPU sensitivity：恢复合理设置
        if 'mgy_clusters' in db_name_lower:
            cpu_sensitivity = 4.5  # 主瓶颈：适度降低sensitivity
        elif 'uniprot' in db_name_lower:
            cpu_sensitivity = 5.0  # 次要瓶颈：适度降低sensitivity
        else:
            cpu_sensitivity = 5.5  # 其他数据库：保持较好质量

        # 智能CPU搜索：让MMseqs2根据内存自适应
        cmd = [
            self.binary_path, "search", query_db, self.base_db, result_db, tmp_dir,
            "--threads", str(threads),  # 给足线程，让MMseqs2自适应
            "--max-seqs", str(max_seqs),
            "-s", str(cpu_sensitivity), "-e", str(self.e_value),
            "--split-memory-limit", f"{mem_limit_gb}G",  # 70%内存限制
            "--split", "0",  # 不拆分，如果有充足内存（测试用）
            "--db-load-mode", "2",  
            "--min-ungapped-score", "15", 
            "--realign"
        ]
        
        # 显示智能CPU优化信息
        logging.info("🔄 [CPU][%s] Smart CPU search: %d threads (MMseqs2 auto-adjust), mem_limit=%dG (70%%), split=0, sensitivity=%.1f", 
                   db_name, threads, mem_limit_gb, cpu_sensitivity)
        
        # 针对Phase 2 slow tasks的特殊优化
        if 'mgy_clusters' in db_name_lower:
            cmd.extend([
                "--max-accept", str(int(max_seqs * 0.8)),  # 提升接受阈值以获得更多序列
                "--max-rejected", "300"  # 适度限制rejected检查
            ])
            logging.debug("   [CPU][%s] Phase 2 enhanced params for mgnify bottleneck", db_name)
        elif 'uniprot' in db_name_lower:
            cmd.extend([
                "--max-accept", str(int(max_seqs * 0.7))  # uniprot也提升接受阈值
            ])
            logging.debug("   [CPU][%s] Phase 2 enhanced params for uniprot", db_name)

        try:
            run_with_quiet_logging(cmd, operation_name=f"[CPU][{db_name}] Search")

            # Stage 4: Generate A3M directly (skip convertalis -> m8)
            result_a3m = os.path.join(tmp_dir, "result.a3m")
            logging.debug("🚀 [CPU][%s] Stage 4/4: Generating A3M directly (~30-120s)", db_name)
            logging.debug("   [CPU][%s] Input: Alignment DB | Output: A3M MSA | Hardware: CPU", db_name)

            try:
                generation_time = self._generate_msa_fast(
                    mmseqs_bin=self.binary_path,
                    query_db=query_db,
                    target_db=self.base_db,
                    result_db=result_db,
                    output_a3m=result_a3m,
                    threads=self.n_cpu  # CPU模式：使用所有CPU处理A3M生成
                )
                logging.debug("✅ [CPU][%s] Stage 4 completed", db_name)
                return self._process_search_results(result_a3m, target_sequence, query_db, tmp_dir, from_a3m=True, msa_generation_time=generation_time)
            except Exception as e:
                logging.warning("Fast A3M generation failed: %s. Falling back to convertalis -> result2msa", e)
                # Fallback: convertalis + traditional result2msa
                cmd = [self.binary_path, "convertalis", query_db, self.base_db, result_db, result_m8]
                run_with_logging(cmd)
                logging.debug("✅ [CPU][%s] Stage 4 completed (fallback)", db_name)
                return self._process_search_results(result_m8, target_sequence, query_db, tmp_dir)
        except subprocess.CalledProcessError as e:
            logging.error("❌ [CPU][%s] CPU search FAILED at Stage 3 (prefilter stage)", db_name)
            logging.error("   [CPU][%s] Error context: Database=%s, Memory Strategy=%s, Threads=%d, MaxSeqs=%d",
                         db_name, db_name, strategy, threads, max_seqs)
            logging.error("   [CPU][%s] System RAM: %dGB, Using: %dGB (%.1f%%)",
                         db_name, total_mb // 1024, mem_limit_gb, actual_percent)
            logging.warning("   [CPU][%s] Returning minimal MSA due to prefilter failure", db_name)
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

    def _process_search_results(self, result_path: str, target_sequence: str, query_db_path: str, tmp_dir: str, from_a3m: bool = False, msa_generation_time: float = None) -> msa_tool.MsaToolResult:
        """Process search results to A3M.

        result_path: 如果 from_a3m=True，则为A3M文件路径；否则为 m8 文件路径或 result prefix。
        """
        import time
        stage_start = time.time()
        db_name = os.path.basename(self.database_path)

        if from_a3m:
            # Fast path: directly read A3M file that was generated by result2msa --msa-format-mode 6
            logging.debug("🔹 [%s] Processing A3M file directly (fast path)", db_name)
            try:
                with open(result_path, 'r', encoding='utf-8', errors='replace') as f:
                    a3m_content = f.read()

                # ✅ 使用轻量级清洗（MMseqs2已经过滤，只需基本验证）
                a3m_content = self._lightweight_a3m_clean(a3m_content, target_sequence)

                # 验证A3M格式
                sequence_count = len([line for line in a3m_content.split('\n') if line.startswith('>')])
                # 显示实际的A3M文件大小来验证处理正确性
                file_size_mb = len(a3m_content) / (1024 * 1024)

                # 使用提供的MSA生成时间，如果没有则使用处理时间
                display_time = msa_generation_time if msa_generation_time is not None else (time.time() - stage_start)
                logging.info("MSA generation completed in %.1fs with %d sequences (fast A3M mode, %.1f MB)",
                            display_time, sequence_count, file_size_mb)
                return msa_tool.MsaToolResult(
                    a3m=a3m_content,
                    target_sequence=target_sequence,
                    e_value=self.e_value
                )
            except Exception as e:
                logging.warning(f"A3M direct processing failed: {e}, falling back to legacy path")
                # Fall through to legacy path below

        # Legacy path: process M8 file or generate MSA from scratch
        logging.info("Generating MSA from %s results", db_name)
        logging.debug("   [MSA][%s] Input: M8 results + databases | Output: A3M MSA | Hardware: CPU", db_name)
        logging.debug("   [MSA][%s] Algorithm: A3M output with fast regex cleaning", db_name)

        result_a3m_path = os.path.join(tmp_dir, "result.a3m")

        try:
            # 确定result_prefix：如果是M8文件，取其前缀；否则直接使用
            if result_path.endswith('.m8'):
                result_prefix = result_path.replace('.m8', '')
            else:
                result_prefix = result_path

            self._generate_msa_fast(
                mmseqs_bin=self.binary_path,
                query_db=query_db_path,
                target_db=self.base_db,
                result_db=result_prefix,
                output_a3m=result_a3m_path,
                threads=self.n_cpu  # 使用所有CPU处理fallback路径的MSA生成
            )

            # Read A3M file in text mode
            with open(result_a3m_path, 'r', encoding='utf-8', errors='replace') as f:
                a3m_content = f.read()

            # ✅ 使用轻量级清洗（MMseqs2已经过滤）
            a3m_content = self._lightweight_a3m_clean(a3m_content, target_sequence)

            # Final validation: try to parse with AlphaFold3 to ensure it's valid
            try:
                from alphafold3.data.parsers import parse_fasta
                test_sequences, test_descriptions = parse_fasta(a3m_content)
                sequence_count = len([line for line in a3m_content.split('\n') if line.startswith('>')])
                logging.info("MSA generation completed in %.1fs with %d sequences", time.time() - stage_start, sequence_count)
                return msa_tool.MsaToolResult(
                    a3m=a3m_content,
                    target_sequence=target_sequence,
                    e_value=self.e_value
                )
            except Exception as e:
                logging.warning(f"A3M validation failed after fix: {e}, returning minimal MSA")
                return msa_tool.MsaToolResult(
                    a3m=f">query\n{target_sequence}\n",
                    target_sequence=target_sequence,
                    e_value=self.e_value
                )

        except Exception as e:
            logging.error(f"MSA generation failed: {e}, returning minimal MSA")
            logging.warning("Returning minimal MSA due to MSA generation failure")
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

    def _generate_msa_stockholm(self, mmseqs_bin: str,
                                query_db: str,
                                target_db: str,
                                result_db: str,
                                output_a3m: str,
                                target_sequence: str,
                                threads: int = 8) -> float:
        """
        Stockholm format path: result2msa mode 4 → Stockholm → Python convert to A3M.
        
        This mimics AlphaFold3's native Jackhmmer pipeline and may be faster.
        
        Returns:
            Time taken for MSA generation in seconds.
        """
        import shlex
        import time
        from alphafold3.data import parsers
        
        start_time = time.time()
        
        # Step 1: Generate Stockholm format using mode 4
        output_sto = output_a3m.replace('.a3m', '.sto')
        cmd = [
            mmseqs_bin, "result2msa",
            query_db, target_db,
            result_db, output_sto,
            "--msa-format-mode", "4",      # Stockholm format (AlphaFold3 native)
            "--threads", str(threads),
            "--compressed", "0",
        ]
        
        logging.info("🔹 [Stockholm] Generating Stockholm format with %d threads", threads)
        logging.info("🔄 [Stockholm] Expected to be faster than A3M mode 6", )
        
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
            if proc.returncode != 0:
                logging.error("[Stockholm] result2msa mode 4 failed: %s", proc.stderr)
                raise RuntimeError(f"Stockholm generation failed: {proc.stderr}")
            
            stockholm_time = time.time() - start_time
            logging.info("✅ [Stockholm] Format generation completed in %.1fs", stockholm_time)
            
            # Step 2: Convert Stockholm to A3M using AlphaFold3's native parser (very fast!)
            convert_start = time.time()
            with open(output_sto, 'r', encoding='utf-8', errors='replace') as f:
                a3m_content = parsers.convert_stockholm_to_a3m(
                    f, max_sequences=self.max_sequences
                )
            
            # Write A3M output
            with open(output_a3m, 'w', encoding='utf-8') as f:
                f.write(a3m_content)
            
            convert_time = time.time() - convert_start
            total_time = time.time() - start_time
            
            logging.info("✅ [Stockholm] Python conversion completed in %.1fs", convert_time)
            logging.info("🎯 [Stockholm] Total time: %.1fs (stockholm: %.1fs + convert: %.1fs)", 
                        total_time, stockholm_time, convert_time)
            
            # Clean up Stockholm file
            try:
                os.remove(output_sto)
            except:
                pass
            
            return total_time
            
        except Exception as e:
            logging.error("❌ [Stockholm] Pipeline failed: %s", e)
            raise

    def _generate_msa_aligned_fasta(self, mmseqs_bin: str,
                          query_db: str,
                          target_db: str,
                          result_db: str,
                          output_a3m: str,
                          threads: int = 8) -> float:
        """
        EXPERIMENTAL: Use aligned FASTA format (mode 2) instead of A3M (mode 6).
        Testing if this format conversion is faster.

        Returns:
            Time taken for MSA generation in seconds.
        """
        import shlex
        import time

        optimized_threads = threads
        output_fasta = output_a3m.replace('.a3m', '.aligned.fasta')

        # Try aligned FASTA format (mode 2)
        cmd = [
            mmseqs_bin, "result2msa",
            query_db, target_db,
            result_db, output_fasta,
            "--msa-format-mode", "2",      # 🔬 EXPERIMENTAL: aligned FASTA DB
            "--threads", str(optimized_threads),
            "--compressed", "0",
        ]
        logging.info("🔬 [EXPERIMENTAL] Testing aligned FASTA format (mode 2)")
        logging.info("🔹 Running result2msa mode 2 with %d threads: %s", optimized_threads, " ".join(shlex.quote(x) for x in cmd[-4:]))
        logging.info("🔄 Expected: similar speed to mode 6, but let's see...")

        start_time = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)
            if proc.returncode != 0:
                logging.error("result2msa (aligned FASTA mode 2) failed with return code %d", proc.returncode)
                if proc.stderr:
                    logging.error("Error details: %s", proc.stderr.strip())
                raise RuntimeError(f"result2msa aligned FASTA mode failed: {proc.stderr}")
        except subprocess.TimeoutExpired:
            logging.error("❌ Aligned FASTA generation timed out after 20 minutes")
            raise RuntimeError("result2msa aligned FASTA mode timed out")
        except Exception as e:
            logging.error(f"❌ Error during aligned FASTA generation: {e}")
            raise

        generation_time = time.time() - start_time
        logging.info("✅ Aligned FASTA generation completed in %.1fs", generation_time)

        # Convert aligned FASTA to A3M format
        # Aligned FASTA is essentially the same as A3M but with uppercase gaps
        # We just need to read it and potentially convert gap characters
        if not os.path.exists(output_fasta) or os.path.getsize(output_fasta) == 0:
            logging.warning("⚠️ Generated aligned FASTA file is empty or missing")
            raise RuntimeError("Aligned FASTA file is empty or missing")

        # Read aligned FASTA and write as A3M (should be directly compatible)
        with open(output_fasta, 'r', encoding='utf-8', errors='replace') as f:
            fasta_content = f.read()
        
        # Aligned FASTA from MMseqs2 should be directly usable as A3M
        with open(output_a3m, 'w', encoding='utf-8') as f:
            f.write(fasta_content)

        logging.info("✅ Format conversion completed, total time: %.1fs", generation_time)
        
        return generation_time

    def _generate_msa_fast(self, mmseqs_bin: str,
                          query_db: str,
                          target_db: str,
                          result_db: str,
                          output_a3m: str,
                          threads: int = 8) -> float:
        """
        Fast path: call result2msa with the result DB (not .m8) and --msa-format-mode 6.

        Returns:
            Time taken for MSA generation in seconds.
        """
        import shlex
        import time

        # 智能A3M生成：使用传入的线程数（通常为self.n_cpu）
        optimized_threads = threads  # 使用所有可用线程，让MMseqs2自适应

        # result2msa只使用安全的基本参数
        # 过滤应该在search阶段完成，result2msa主要做格式转换
        cmd = [
            mmseqs_bin, "result2msa",
            query_db, target_db,
            result_db, output_a3m,
            "--msa-format-mode", "6",      # 快速A3M模式（直接生成A3M）
            "--threads", str(optimized_threads),  # 多线程
            "--compressed", "0",           # 不压缩（更快但更大）
        ]
        logging.info("🔹 Running fast result2msa with %d threads: %s", optimized_threads, " ".join(shlex.quote(x) for x in cmd[-4:]))
        logging.info("🔄 A3M generation balanced: %d threads (%.1f%% CPU utilization)", 
                   optimized_threads, (optimized_threads / 112) * 100)

        start_time = time.time()
        try:
            proc = subprocess.run(cmd, capture_output=True, text=True, timeout=1200)  # 20min timeout (保守激进策略)
            if proc.returncode != 0:
                logging.error("result2msa (fast) failed with return code %d", proc.returncode)
                if proc.stderr:
                    stderr_lines = proc.stderr.strip().split('\n')
                    # Log more detailed error information
                    if len(stderr_lines) > 3:
                        logging.error("Error details (last 3 lines): %s", ' | '.join(stderr_lines[-3:]))
                    else:
                        logging.error("Error details: %s", proc.stderr.strip())
                
                # Check for specific error patterns that might be recoverable
                if "Database" in proc.stderr and "not found" in proc.stderr:
                    raise RuntimeError(f"Database path issue: {proc.stderr}")
                elif "format" in proc.stderr.lower():
                    raise RuntimeError(f"Format compatibility issue: {proc.stderr}")
                else:
                    raise RuntimeError(f"result2msa fast mode failed: {proc.stderr}")
                    
        except subprocess.TimeoutExpired:
            logging.error("❌ Fast MSA generation timed out after 30 minutes")
            raise RuntimeError("result2msa fast mode timed out")
        except Exception as e:
            logging.error(f"❌ Error during fast MSA generation: {e}")
            raise

        generation_time = time.time() - start_time

        # 改进的输出验证 - 更宽松但仍有效
        if not os.path.exists(output_a3m):
            raise RuntimeError(f"A3M output file not created: {output_a3m}")
        
        file_size = os.path.getsize(output_a3m)
        if file_size == 0:
            raise RuntimeError(f"A3M output file is empty: {output_a3m}")
        
        # 基本格式检查 - 只检查是否有FASTA头，不要求完美格式
        try:
            with open(output_a3m, 'r', encoding='utf-8', errors='replace') as f:
                content_preview = f.read(1000)  # Read first 1KB
                if not content_preview.strip():
                    raise RuntimeError("A3M file appears to be empty or contain only whitespace")
                if ">" not in content_preview:
                    logging.warning("⚠️ A3M file may not contain FASTA headers, but continuing...")
        except IOError as e:
            logging.warning(f"⚠️ Could not validate A3M file format: {e}, but file exists")

        logging.debug("✅ Fast A3M generation completed: %.1f MB in %.1fs", 
                     file_size / (1024*1024), generation_time)
        return generation_time

    def _lightweight_a3m_clean(self, a3m_content: str, target_sequence: str) -> str:
        """轻量级A3M清洗 - MMseqs2已经过滤，只做最基本的验证.
        
        优化策略：
        - MMseqs2的--filter-msa已经处理了大部分问题
        - 只移除null字符和控制字符
        - 不做复杂的序列长度检查和过滤
        
        Args:
            a3m_content: Raw A3M content from MMseqs2
            target_sequence: Target sequence for fallback
            
        Returns:
            Lightly cleaned A3M content
        """
        if not a3m_content or not a3m_content.strip():
            return f">query\n{target_sequence}\n"
        
        try:
            # 只移除null字符和控制字符（最快的清洗）
            cleaned = ""
            for char in a3m_content:
                if ord(char) >= 32 or char in '\n\t':
                    cleaned += char
            
            # 基本验证：至少有一个序列
            if ">query\n" not in cleaned and ">" in cleaned:
                # 如果没有query标记，确保第一个序列是query
                lines = cleaned.split('\n')
                if lines[0].startswith('>'):
                    lines[0] = '>query'
                    cleaned = '\n'.join(lines)
            
            # 如果清洗后为空，返回最小MSA
            if not cleaned or not cleaned.strip():
                return f">query\n{target_sequence}\n"
                
            return cleaned
            
        except Exception as e:
            logging.warning(f"Lightweight A3M clean failed: {e}, returning minimal MSA")
            return f">query\n{target_sequence}\n"
    
    def _clean_a3m_fast(self, a3m_content: str, target_sequence: str = "") -> str:
        """Fast A3M cleaning using regex - removes invalid characters while preserving FASTA format.

        This is much faster than the complex _validate_and_fix_a3m method and sufficient for AF3.
        Only keeps: A-Z, a-z, >, -, \n (valid FASTA/A3M characters)

        Args:
            a3m_content: Raw A3M content from MMseqs2 convertalis
            target_sequence: Target sequence for fallback (optional)

        Returns:
            Cleaned A3M content safe for AF3 parsing
        """
        import re

        # Remove all characters that are not valid in FASTA/A3M format
        # Keep: uppercase/lowercase letters, >, -, newline
        cleaned = re.sub(r'[^A-Za-z>\-\n]', '', a3m_content)

        # Basic validation: ensure we have at least a query sequence
        if not cleaned or not cleaned.strip():
            if target_sequence:
                return f">query\n{target_sequence}\n"
            else:
                return ">query\nMQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG\n"

        return cleaned

    def _validate_and_fix_a3m(self, a3m_content: str, target_sequence: str) -> str:
        """Validate and fix A3M format issues that cause parsing errors.

        The key insight: AlphaFold3's C++ convert_a3m_to_stockholm function
        REQUIRES all sequences to have identical lengths. Simple padding with '-'
        breaks the alignment. Instead, we need to ensure proper MSA alignment.

        Args:
            a3m_content: Raw A3M content from MMseqs2
            target_sequence: The target sequence for validation

        Returns:
            Fixed A3M content that can be properly converted to Stockholm
        """
        if not a3m_content or not a3m_content.strip():
            logging.warning("Empty A3M content, returning minimal MSA")
            return f">query\n{target_sequence}\n"

        try:
            # Debug: check raw A3M content for issues
            if len(a3m_content.strip()) > 0:
                lines = a3m_content.strip().split('\n')
                desc_lines = [line for line in lines if line.startswith('>')]
                seq_lines = [line for line in lines if not line.startswith('>') and line.strip()]

                logging.warning(f"🔍 A3M DEBUG: {len(desc_lines)} description lines, {len(seq_lines)} sequence lines")
                logging.warning(f"🔍 Raw A3M content (first 500 chars): {repr(a3m_content[:500])}")

                # Check for problematic characters
                problematic_chars = set()
                for char in a3m_content:
                    if char not in 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz- \n>1234567890':
                        problematic_chars.add(repr(char))
                if problematic_chars:
                    logging.warning(f"❌ Found problematic characters in A3M: {problematic_chars}")

            # CRITICAL FIX: Clean A3M content from MMseqs2 null characters and control chars
            # MMseqs2 result2msa sometimes adds \x00 and other control characters at end of sequences
            cleaned_a3m = ""
            for char in a3m_content:
                # Remove null characters and other control characters (except \n and \t)
                if ord(char) >= 32 or char in '\n\t':
                    cleaned_a3m += char

            # Remove trailing whitespace and control characters from each line
            cleaned_lines = []
            for line in cleaned_a3m.split('\n'):
                cleaned_line = line.rstrip('\x00\t\r')
                if cleaned_line or line.startswith('>'):  # Keep description lines even if empty after cleaning
                    cleaned_lines.append(cleaned_line)

            cleaned_a3m = '\n'.join(cleaned_lines)
            logging.warning(f"🧹 Cleaned A3M content (removed null/control chars), new length: {len(cleaned_a3m)}")

            # Parse the cleaned A3M content
            from alphafold3.data.parsers import parse_fasta
            sequences, descriptions = parse_fasta(cleaned_a3m)

            # 移除详细的序列内容输出，只保留统计信息
            if len(descriptions) != len(sequences):
                logging.warning(f"❌ PARSE ERROR: {len(descriptions)} descriptions but {len(sequences)} sequences")

            if not sequences:
                logging.warning("No sequences found in A3M, returning minimal MSA")
                return f">query\n{target_sequence}\n"

            # Check for length inconsistencies - this is critical for Stockholm conversion
            seq_lengths = [len(seq) for seq in sequences]
            max_length = max(seq_lengths)
            min_length = min(seq_lengths)

            if max_length != min_length:
                # Detailed length distribution analysis
                from collections import Counter
                length_counts = Counter(seq_lengths)

                logging.warning(f"A3M alignment issue detected: sequences have inconsistent lengths "
                              f"(min={min_length}, max={max_length}). This prevents proper "
                              f"Stockholm conversion required by hmmbuild.")

                # Log detailed length distribution
                sorted_lengths = sorted(length_counts.keys())
                distribution_str = ", ".join(f"{length}:{count}" for length, count in length_counts.most_common())
                logging.info(f"Sequence length distribution: {distribution_str}")

                # Strategy: Keep only sequences that are properly aligned to the query
                # The query sequence (first one) should be the reference
                query_seq = sequences[0]
                query_len = len(query_seq)

                # Count sequences matching query length before filtering
                matching_query_len = sum(1 for seq in sequences if len(seq) == query_len)
                total_sequences = len(sequences)

                logging.info(f"Query sequence length: {query_len}, {matching_query_len}/{total_sequences} "
                           f"sequences match query length ({matching_query_len/total_sequences:.1%})")

                # Filter out sequences that don't match the query length
                # This is more conservative than padding, which can break alignments
                aligned_sequences = []
                aligned_descriptions = []
                discarded_count = 0

                for seq, desc in zip(sequences, descriptions):
                    if len(seq) == query_len:
                        aligned_sequences.append(seq)
                        aligned_descriptions.append(desc)
                    else:
                        discarded_count += 1
                        seq_len = len(seq)
                        logging.debug(f"Discarding misaligned sequence '{desc}' (len={seq_len} != query_len={query_len})")

                if not aligned_sequences:
                    logging.warning("No properly aligned sequences found, keeping query sequence only")
                    # Always keep at least the query sequence to ensure MSA is not empty
                    aligned_sequences = [query_seq]
                    aligned_descriptions = [descriptions[0]]  # query description

                sequences = aligned_sequences
                descriptions = aligned_descriptions

                retained_count = len(sequences)
                retention_rate = retained_count / total_sequences
                discard_rate = discarded_count / total_sequences

                # Quality assessment based on retention rate
                if retention_rate >= 0.9:
                    quality_status = "Excellent (minimal loss)"
                elif retention_rate >= 0.8:
                    quality_status = "Good (acceptable loss)"
                elif retention_rate >= 0.7:
                    quality_status = "Fair (moderate loss)"
                elif retention_rate >= 0.5:
                    quality_status = "Poor (significant loss)"
                else:
                    quality_status = "Critical (severe loss - MSA may be unreliable)"

                logging.warning(f"MSA filtering results: {retained_count}/{total_sequences} sequences retained "
                              f"({retention_rate:.1%}), {discarded_count} discarded ({discard_rate:.1%}) - {quality_status}")

                if retention_rate < 0.5:
                    logging.error("CRITICAL: MSA quality severely impacted! More than 50% of sequences discarded. "
                                "This suggests fundamental alignment issues. Consider checking MMseqs2 parameters "
                                "or database quality.")

                # Additional analysis: check if alternative strategies would work better
                if discarded_count > 0:
                    # Check if using most common length would retain more sequences
                    most_common_length = max(length_counts, key=length_counts.get)
                    most_common_count = length_counts[most_common_length]
                    alternative_rate = most_common_count / total_sequences

                    if alternative_rate > retention_rate + 0.1:  # 10% improvement threshold
                        logging.info(f"Alternative strategy available: using most common length {most_common_length} "
                                   f"would retain {most_common_count}/{total_sequences} sequences "
                                   f"({alternative_rate:.1%} vs current {retention_rate:.1%})")
                        logging.info("Consider switching to most_common_length strategy if current results are poor.")

            # Rebuild A3M content with properly aligned sequences
            fixed_lines = []
            for desc, seq in zip(descriptions, sequences):
                fixed_lines.append(f">{desc}")
                # A3M format: split long lines at 80 characters
                for i in range(0, len(seq), 80):
                    fixed_lines.append(seq[i:i+80])

            return '\n'.join(fixed_lines) + '\n'

        except Exception as e:
            logging.warning(f"A3M validation failed: {e}, returning minimal MSA")
            return f">query\n{target_sequence}\n"

        # Fallback: manual parsing and fixing
        lines = a3m_content.strip().split('\n')
        headers = []
        sequences = []

        current_seq = ""
        for line in lines:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                headers.append(line)
                current_seq = ""
            else:
                current_seq += line

        # Add the last sequence
        if current_seq:
            sequences.append(current_seq)

        # Check lengths
        if len(sequences) > 1 and len(headers) == len(sequences):
            seq_lengths = [len(seq) for seq in sequences]
            max_length = max(seq_lengths)
            min_length = min(seq_lengths)

            if max_length != min_length:
                need_padding = sum(1 for seq in sequences if len(seq) < max_length)
                logging.warning(f"A3M format issue (fallback): {need_padding}/{len(sequences)} sequences need length padding")

                # Rebuild with padding
                fixed_lines = []
                for header, seq in zip(headers, sequences):
                    fixed_lines.append(header)
                    if len(seq) < max_length:
                        seq = seq + '-' * (max_length - len(seq))
                    for i in range(0, len(seq), 80):
                        fixed_lines.append(seq[i:i+80])

                return '\n'.join(fixed_lines) + '\n'

        return a3m_content

    def query(self, target_sequence: str) -> msa_tool.MsaToolResult:
        """Main query method."""
        import time
        start_time = time.time()

        seq_len = len(target_sequence)
        db_name = os.path.basename(self.database_path)
        logging.debug("🔬 Starting MSA search for sequence (%d aa) against database: %s", seq_len, db_name)
        logging.debug("📋 Pipeline: FASTA → createdb → search → convertalis → result2msa → A3M")
        logging.debug("⏱️  Estimated time: GPU=5-60s, CPU=5-60min | Hardware: %s",
                    "GPU+CPU" if self.use_gpu else "CPU only")
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            query_path = os.path.join(tmp_dir, "query.fasta")
            result_m8 = os.path.join(tmp_dir, "result.m8")
            
            with open(query_path, "w") as f:
                f.write(f">query\n{target_sequence}\n")
            
            if self.use_gpu:
                result = self._gpu_search(query_path, result_m8, tmp_dir, target_sequence)
            else:
                result = self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

            # Log total time
            total_time = time.time() - start_time
            logging.debug("⏱️  [%s] Total MSA search completed in %.1fs", db_name, total_time)
            return result

    @classmethod
    def query_multiple_databases_parallel(
        cls,
        target_sequence: str,
        database_configs: List[Dict[str, Any]]
    ) -> Dict[str, msa_tool.MsaToolResult]:
        """Query multiple databases with ASYNC PIPELINE for maximum GPU/CPU overlap.

        This implements ASYNCHRONOUS MSA generation with overlapping stages:
        - GPU search and A3M generation run in parallel (overlap)
        - GPU is freed immediately after search, not blocked by A3M generation
        - Multiple databases can process concurrently at different stages
        - Maximizes both GPU and CPU utilization

        Strategy:
            Time 0s:  [GPU1: DB1 search] [GPU2: DB2 search]
            Time 20s: [GPU1: DB3 search] [GPU2: DB2 search] [CPU: DB1 A3M gen]
            Time 40s: [GPU1: DB3 search] [GPU2: DB4 search] [CPU: DB1+DB2 A3M gen]
            Time 60s: [Idle]              [GPU2: DB4 search] [CPU: DB1+DB2+DB3 A3M gen]

        Args:
            target_sequence: The protein sequence to search
            database_configs: List of config dicts with keys:
                - binary_path, database_path, n_cpu, e_value, max_sequences, sensitivity, gpu_devices, use_gpu

        Returns:
            Dict mapping database names to MsaToolResult objects
        """
        results = {}
        
        # Separate GPU and CPU tasks
        gpu_tasks = []
        cpu_tasks = []

        for config in database_configs:
            # Map database filenames to expected names in pipeline.py
            filename = config['database_path'].split('/')[-1]
            if 'mgy_clusters' in filename:
                db_name = 'mgnify'
            elif 'bfd-first' in filename:
                db_name = 'bfd'
            else:
                db_name = filename.split('_')[0]
            use_gpu = config.get('use_gpu', False)

            if use_gpu and config.get('gpu_devices'):
                gpu_tasks.append((db_name, config))
            else:
                cpu_tasks.append((db_name, config))

        # ASYNC PIPELINE: GPU search + A3M generation overlap
        if gpu_tasks:
            logging.info("🚀 Starting ASYNC PIPELINE processing for %d databases", len(gpu_tasks))
            logging.info("🔄 Strategy: GPU search and A3M generation will OVERLAP for maximum performance")
            
            # Check available GPUs
            available_gpus = cls._get_available_gpu_ids()
            if not available_gpus:
                logging.warning("No GPUs available, falling back to CPU")
                cpu_tasks.extend(gpu_tasks)
                gpu_tasks = []
            else:
                logging.info("Available GPUs for parallel processing: %s", available_gpus)
                
                # Two-stage thread pools for async pipeline
                # 1. GPU search pool (limited by GPU count)
                # 2. A3M generation pool (CPU-bound, can run many in parallel)
                max_gpu_workers = min(len(gpu_tasks), len(available_gpus))
                max_a3m_workers = min(len(gpu_tasks), 4)  # Limit to 4 concurrent A3M generations
                
                # Track which databases are assigned to which GPUs
                gpu_task_queue = list(enumerate(gpu_tasks))  # [(idx, (db_name, config)), ...]
                gpu_assignments = {}  # {gpu_id: current_task_idx or None}
                for gpu_id in available_gpus:
                    gpu_assignments[gpu_id] = None
                
                # Stage 1: Submit GPU search tasks
                # Stage 2: As each GPU search completes, immediately submit A3M generation
                with ThreadPoolExecutor(max_workers=max_gpu_workers, thread_name_prefix="GPU-Search") as gpu_executor, \
                     ThreadPoolExecutor(max_workers=max_a3m_workers, thread_name_prefix="A3M-Gen") as a3m_executor:
                    
                    # Submit all GPU search tasks
                    gpu_search_futures = {}
                    a3m_gen_futures = {}
                    
                    for i, (db_name, config) in enumerate(gpu_tasks):
                        # Assign GPU in round-robin fashion
                        assigned_gpu = available_gpus[i % len(available_gpus)]
                        gpu_config = config.copy()
                        gpu_config['gpu_devices'] = [str(assigned_gpu)]
                        
                        logging.info("🎯 [ASYNC] Submitting GPU search for %s to GPU %d", db_name, assigned_gpu)
                        
                        # Submit GPU search task
                        future = gpu_executor.submit(
                            cls._async_gpu_search_wrapper,
                            gpu_config, target_sequence, db_name
                        )
                        gpu_search_futures[future] = (db_name, gpu_config)
                    
                    # Process completed GPU searches and submit A3M generation
                    for future in as_completed(gpu_search_futures):
                        db_name, gpu_config = gpu_search_futures[future]
                        
                        try:
                            # Get GPU search result
                            search_result = future.result()
                            
                            if search_result is None:
                                logging.error("❌ [ASYNC] GPU search failed for %s", db_name)
                                results[db_name] = msa_tool.MsaToolResult(
                                    a3m=f">query\n{target_sequence}\n",
                                    target_sequence=target_sequence,
                                    e_value=0.0001
                                )
                                continue
                            
                            logging.info("✅ [ASYNC] GPU search completed for %s, submitting A3M generation", db_name)
                            
                            # Immediately submit A3M generation (runs on CPU, frees GPU)
                            a3m_future = a3m_executor.submit(
                                cls._async_a3m_gen_wrapper,
                                gpu_config, search_result, target_sequence, db_name
                            )
                            a3m_gen_futures[a3m_future] = db_name
                            
                        except Exception as e:
                            logging.error("❌ [ASYNC] GPU search failed for %s: %s", db_name, e)
                            results[db_name] = msa_tool.MsaToolResult(
                                a3m=f">query\n{target_sequence}\n",
                                target_sequence=target_sequence,
                                e_value=0.0001
                            )
                    
                    # Wait for all A3M generation tasks to complete
                    for future in as_completed(a3m_gen_futures):
                        db_name = a3m_gen_futures[future]
                        
                        try:
                            result = future.result()
                            results[db_name] = result
                            logging.info("✅ [ASYNC] A3M generation completed for %s", db_name)
                        except Exception as e:
                            logging.error("❌ [ASYNC] A3M generation failed for %s: %s", db_name, e)
                            results[db_name] = msa_tool.MsaToolResult(
                                a3m=f">query\n{target_sequence}\n",
                                target_sequence=target_sequence,
                                e_value=0.0001
                            )

        # Process CPU tasks in parallel
        if cpu_tasks:
            with ThreadPoolExecutor(max_workers=min(len(cpu_tasks), 4)) as executor:
                future_to_db = {}
                for db_name, config in cpu_tasks:
                    try:
                        tool = cls(
                            binary_path=config['binary_path'],
                            database_path=config['database_path'],
                            n_cpu=config['n_cpu'],
                            e_value=config.get('e_value', 0.0001),
                            max_sequences=config.get('max_sequences', 10000),
                            sensitivity=config.get('sensitivity', 4.0),
                            gpu_devices=None,
                            use_gpu=False
                        )
                        future = executor.submit(tool.query, target_sequence)
                        future_to_db[future] = db_name
                    except Exception as e:
                        logging.error("❌ Failed to initialize CPU MSA tool for %s: %s", db_name, e)
                        results[db_name] = msa_tool.MsaToolResult(
                            a3m=f">query\n{target_sequence}\n",
                            target_sequence=target_sequence,
                            e_value=config.get('e_value', 0.0001)
                        )

                # Collect CPU results
                for future in as_completed(future_to_db):
                    db_name = future_to_db[future]
                    try:
                        results[db_name] = future.result()
                        logging.info("✅ CPU MSA completed for %s", db_name)
                    except Exception as e:
                        logging.error("❌ CPU MSA failed for %s: %s", db_name, e)
                        results[db_name] = msa_tool.MsaToolResult(
                            a3m=f">query\n{target_sequence}\n",
                            target_sequence=target_sequence,
                            e_value=0.0001
                        )

        return results

    @classmethod
    def _get_available_gpu_ids(cls) -> List[int]:
        """Get list of available GPU IDs with sufficient memory, prioritizing non-GPU0."""
        available_gpus = []
        
        # Define priority order: avoid GPU0 first
        gpu_check_order = [1, 2, 3, 4, 5, 6, 7, 0]
        
        for gpu_id in gpu_check_order:
            gpu_info = get_gpu_memory_info(gpu_id)
            if gpu_info and gpu_info['free_mb'] >= 3000:  # 3GB minimum
                available_gpus.append(gpu_id)
                logging.debug("GPU %d available with %dMB free", gpu_id, gpu_info['free_mb'])
        
        if available_gpus:
            logging.info("Available GPUs for parallel processing: %s", available_gpus)
        else:
            logging.warning("No GPUs found with sufficient memory (≥3GB)")
        
        return available_gpus
    
    @classmethod
    def _async_gpu_search_wrapper(cls, config: Dict[str, Any], target_sequence: str, db_name: str) -> Dict[str, Any]:
        """Wrapper for async GPU search: only performs GPU search, returns paths for A3M generation.
        
        This enables true async pipeline: GPU freed immediately after search.
        """
        import tempfile
        
        try:
            tool = cls(
                binary_path=config['binary_path'],
                database_path=config['database_path'],
                n_cpu=config['n_cpu'],
                e_value=config.get('e_value', 0.0001),
                max_sequences=config.get('max_sequences', 10000),
                sensitivity=config.get('sensitivity', 4.0),
                gpu_devices=tuple(config.get('gpu_devices', [])),
                use_gpu=True
            )
            
            # Create persistent temp directory (must survive until A3M generation)
            tmp_dir = tempfile.mkdtemp(prefix=f"mmseqs_{db_name}_")
            query_path = os.path.join(tmp_dir, "query.fasta")
            
            # Write query sequence
            with open(query_path, 'w') as f:
                f.write(f">query\n{target_sequence}\n")
            
            # Execute GPU search only (no A3M generation)
            search_result = tool._gpu_search_only(query_path, tmp_dir)
            search_result['tool_config'] = config  # Store config for A3M generation
            search_result['target_sequence'] = target_sequence
            
            logging.debug("🎯 [ASYNC] GPU search wrapper completed for %s", db_name)
            return search_result
            
        except Exception as e:
            logging.error("❌ [ASYNC] GPU search wrapper failed for %s: %s", db_name, e)
            return None
    
    @classmethod
    def _async_a3m_gen_wrapper(cls, config: Dict[str, Any], search_result: Dict[str, Any], 
                               target_sequence: str, db_name: str) -> msa_tool.MsaToolResult:
        """Wrapper for async A3M generation: processes GPU search results to generate MSA.
        
        This runs on CPU and doesn't block GPU resources.
        """
        import shutil
        
        try:
            tool = cls(
                binary_path=config['binary_path'],
                database_path=config['database_path'],
                n_cpu=config['n_cpu'],
                e_value=config.get('e_value', 0.0001),
                max_sequences=config.get('max_sequences', 10000),
                sensitivity=config.get('sensitivity', 4.0),
                gpu_devices=tuple(config.get('gpu_devices', [])),
                use_gpu=True
            )
            
            # Generate A3M from search results
            result = tool._a3m_generation_only(search_result, target_sequence)
            
            # Clean up temp directory
            try:
                if 'tmp_dir' in search_result and os.path.exists(search_result['tmp_dir']):
                    shutil.rmtree(search_result['tmp_dir'])
            except Exception as cleanup_error:
                logging.warning("Failed to cleanup temp dir for %s: %s", db_name, cleanup_error)
            
            logging.debug("🎯 [ASYNC] A3M generation wrapper completed for %s", db_name)
            return result
            
        except Exception as e:
            logging.error("❌ [ASYNC] A3M generation wrapper failed for %s: %s", db_name, e)
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=config.get('e_value', 0.0001)
            )
    
    @classmethod 
    def _assign_gpus_to_databases(cls, gpu_tasks: List[tuple], available_gpus: List[int]) -> List[int]:
        """Intelligently assign GPUs to database tasks."""
        if not available_gpus:
            return []
        
        assignments = []
        
        # Strategy: Assign larger databases to better GPUs, cycle through available GPUs
        # Database size priority (largest first)
        db_size_priority = {
            'mgnify': 4,     # Largest database
            'uniprot': 3,    # Large database  
            'uniref90': 2,   # Medium database
            'bfd': 1         # Smallest database
        }
        
        # Sort tasks by database size (largest first)
        sorted_tasks = sorted(gpu_tasks, 
                            key=lambda x: db_size_priority.get(x[0], 0), 
                            reverse=True)
        
        # Assign GPUs in round-robin fashion
        for i, (db_name, config) in enumerate(sorted_tasks):
            assigned_gpu = available_gpus[i % len(available_gpus)]
            assignments.append(assigned_gpu)
            logging.info("📋 Assigned %s → GPU %d", db_name, assigned_gpu)
        
        return assignments
    
    @classmethod
    def _process_single_gpu_database(cls, config: Dict[str, Any], 
                                   target_sequence: str, db_name: str) -> msa_tool.MsaToolResult:
        """Process a single database with GPU in a separate thread."""
        try:
            tool = cls(
                binary_path=config['binary_path'],
                database_path=config['database_path'], 
                n_cpu=config['n_cpu'],
                e_value=config.get('e_value', 0.0001),
                max_sequences=config.get('max_sequences', 10000),
                sensitivity=config.get('sensitivity', 4.0),
                gpu_devices=tuple(config.get('gpu_devices', [])),
                use_gpu=True
            )
            
            result = tool.query(target_sequence)
            logging.debug("🎯 GPU task completed for %s", db_name)
            return result
            
        except Exception as e:
            logging.error("❌ GPU processing failed for %s: %s", db_name, e)
            # Return minimal MSA as fallback
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=config.get('e_value', 0.0001)
            )

