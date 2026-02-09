"""Simple MMseqs2 MSA tool with memory control."""

import os
import subprocess
import tempfile
from typing import Sequence

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
        self.use_gpu = use_gpu and gpu_devices is not None and len(gpu_devices) > 0
        self.gpu_devices = gpu_devices or ()

        # Verify binary
        subprocess_utils.check_binary_exists(path=self.binary_path, name='MMseqs2')

        # Setup database
        if not os.path.exists(self.database_path):
            raise ValueError(f"Database not found: {self.database_path}")

        db_dir = os.path.dirname(self.database_path)
        db_name = os.path.splitext(os.path.basename(self.database_path))[0]
        self.mmseqs_dir = os.path.join(db_dir, f"{db_name}_mmseqs2")
        os.makedirs(self.mmseqs_dir, exist_ok=True)
        self.base_db = os.path.join(self.mmseqs_dir, db_name)

        # Create base database if needed
        if not os.path.exists(self.base_db + ".dbtype"):
            cmd = [self.binary_path, 'createdb', self.database_path, self.base_db, '--dbtype', '1']
            run_with_logging(cmd)

        # GPU setup
        if self.use_gpu:
            self.gpu_index_dir = os.path.join(self.mmseqs_dir, f"{db_name}_gpu_index")
            os.makedirs(self.gpu_index_dir, exist_ok=True)
            self.gpu_db = os.path.join(self.gpu_index_dir, f"{db_name}_gpu")
            
            gpu_files = [self.gpu_db + ext for ext in [".dbtype", ".index", ".lookup", "_h"]]
            if not all(os.path.exists(f) for f in gpu_files):
                cmd = [self.binary_path, 'makepaddedseqdb', self.base_db, self.gpu_db]
                run_with_logging(cmd)
                with tempfile.TemporaryDirectory() as tmp_dir:
                    cmd = [self.binary_path, 'createindex', self.gpu_db, tmp_dir, '--threads', str(self.n_cpu)]
                    run_with_logging(cmd)

    def _gpu_search(self, query_path: str, result_m8: str, tmp_dir: str, target_sequence: str) -> msa_tool.MsaToolResult:
        """GPU search with memory check."""
        import time
        db_name = os.path.basename(self.database_path)
        db_short = db_name.split('_')[0] if '_' in db_name else db_name[:10]
        logging.info("🚀 [%s] Starting GPU search", db_short)
        logging.info("📊 [%s] Database: %s | Size: ~%s | Hardware: GPU (CUDA)",
                    db_short, db_name, _estimate_db_size(db_name))

        if not self.gpu_db:
            logging.warning("GPU database not available, falling back to CPU")
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)
        
        # Stage 1: Create query database (CPU, ~2-10s)
        stage_start = time.time()
        logging.info("📝 [%s] Stage 1/5: Creating query database (CPU, ~2-10s)", db_short)
        logging.info("   [%s] Input: FASTA sequence | Output: MMseqs2 DB | Hardware: CPU", db_short)
        query_db = os.path.join(tmp_dir, "query_db")
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)
        logging.info("✅ [%s] Stage 1 completed in %.1fs", db_short, time.time() - stage_start)
        
        # Stage 2: GPU memory check and selection (~1-5s)
        stage_start = time.time()
        logging.info("🔍 [%s] Stage 2/5: Checking GPU memory and selecting device (~1-5s)", db_short)
        logging.info("   [%s] Hardware: GPU | Requirement: ≥4000MB free memory", db_short)
        gpu_id = None
        for device_id in self.gpu_devices:
            gpu_info = get_gpu_memory_info(int(device_id))
            if gpu_info:
                logging.info("   [%s] GPU %s: %sMB free / %sMB total",
                            db_short, device_id, gpu_info['free_mb'], gpu_info['total_mb'])
                if gpu_info['free_mb'] >= 4000:  # Lower threshold to allow GPU usage
                    gpu_id = device_id
                    break

        if gpu_id is None:
            logging.warning("[%s] No GPU has sufficient memory (need 4000MB+ free), using CPU", db_short)
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

        logging.info("✅ [%s] Selected GPU %s for search", db_short, gpu_id)
        logging.info("✅ [%s] Stage 2 completed in %.1fs", db_short, time.time() - stage_start)

        # Stage 3: GPU search (main bottleneck, ~10-300s)
        stage_start = time.time()
        result_db = os.path.join(tmp_dir, "result")
        env = os.environ.copy()
        env["CUDA_VISIBLE_DEVICES"] = str(gpu_id)

        # GPU memory allocation and max sequences for AF3 small databases - optimized for speed
        db_name = os.path.basename(self.database_path)
        if 'uniprot_all' in db_name.lower():
            gpu_memory_limit = "512M"  # Reduced for speed
            gpu_max_seqs = 100
        elif 'mgy_clusters' in db_name.lower():
            gpu_memory_limit = "256M"  # Reduced for speed
            gpu_max_seqs = 50
        else:
            gpu_memory_limit = "128M"  # Minimal for speed
            gpu_max_seqs = 30

        cmd = [
            self.binary_path, "search", query_db, self.gpu_db, result_db, tmp_dir,
            "--threads", str(min(self.n_cpu, 4)), "--max-seqs", str(gpu_max_seqs), "-s", str(self.sensitivity),
            "-e", str(self.e_value), "--gpu", "1", "--split-memory-limit", gpu_memory_limit,
            "--min-ungapped-score", "10"
        ]

        try:
            run_with_quiet_logging(cmd, env=env, operation_name="GPU search")
            logging.info("✅ [%s] Stage 3 completed in %.1fs", db_short, time.time() - stage_start)

            # Stage 4: Convert results (CPU, ~5-30s)
            stage_start = time.time()
            logging.info("📋 [%s] Stage 4/5: Converting results to M8 format (~5-30s)", db_short)
            logging.info("   [%s] Input: Alignment DB | Output: BLAST-like TSV | Hardware: CPU", db_short)
            cmd = [self.binary_path, "convertalis", query_db, self.gpu_db, result_db, result_m8]
            run_with_logging(cmd)
            logging.info("✅ [%s] Stage 4 completed in %.1fs", db_short, time.time() - stage_start)

            return self._process_search_results(result_m8, target_sequence, query_db, tmp_dir)
        except subprocess.CalledProcessError:
            logging.warning("GPU search failed, falling back to CPU")
            return self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

    def _cpu_search(self, query_path: str, result_m8: str, tmp_dir: str, target_sequence: str) -> msa_tool.MsaToolResult:
        """CPU search with memory limit."""
        db_name = os.path.basename(self.database_path)
        logging.info("🐌 Starting CPU search for %s", db_name)
        logging.info("📊 Database: %s | Size: ~%s | Hardware: CPU only",
                    db_name, _estimate_db_size(db_name))

        # Stage 1: Create query database (CPU, ~2-10s)
        logging.info("📝 Stage 1/4: Creating query database (CPU, ~2-10s)")
        logging.info("   Input: FASTA sequence | Output: MMseqs2 DB | Hardware: CPU")
        query_db = os.path.join(tmp_dir, "query_db")
        cmd = [self.binary_path, "createdb", query_path, query_db]
        run_with_logging(cmd)

        # Stage 2: Memory limit calculation (~1s)
        logging.info("⚙️  Stage 2/4: Calculating memory limits (~1s)")
        db_name = os.path.basename(self.database_path)
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        total_kb = int(line.split()[1])
                        total_mb = total_kb // 1024
                        break
        except:
            total_mb = 256 * 1024  # Assume 256GB if can't read

        # Database-specific memory allocation and max sequences for AF3 small databases - heavily optimized for speed
        if 'uniprot_all' in db_name.lower():
            mem_limit_gb = min(4, max(1, int(total_mb * 0.008) // 1024))  # Up to 4GB for uniprot_all (reduced)
            mem_percent = 0.8
            max_seqs = 100  # Heavily reduced for speed
        elif 'mgy_clusters' in db_name.lower():
            mem_limit_gb = min(2, max(1, int(total_mb * 0.006) // 1024))  # Up to 2GB for mgy_clusters (reduced)
            mem_percent = 0.6
            max_seqs = 50   # Heavily reduced for speed
        else:
            mem_limit_gb = min(2, max(1, int(total_mb * 0.005) // 1024))  # Up to 2GB for small DBs (reduced)
            mem_percent = 0.5
            max_seqs = 30   # Heavily reduced for speed

        logging.info("   System RAM: ~%dGB | Using: %dGB (%d%%) | Algorithm: k-mer + ungapped + gapped",
                    total_mb // 1024, mem_limit_gb, mem_percent)

        # Stage 3: CPU search (main bottleneck, ~5-60min)
        logging.info("🔍 Stage 3/4: Executing CPU homology search (~5-60min)")
        logging.info("   Input: Query DB + Target DB | Output: Alignment DB | Hardware: CPU")
        logging.info("   Algorithm: k-mer matching + ungapped prefilter + gapped alignment")

        result_db = os.path.join(tmp_dir, "result")
        cmd = [
            self.binary_path, "search", query_db, self.base_db, result_db, tmp_dir,
            "--threads", str(min(self.n_cpu, 4)), "--max-seqs", str(max_seqs),
            "-s", str(self.sensitivity), "-e", str(self.e_value),
            "--split-memory-limit", f"{mem_limit_gb}G", "--min-ungapped-score", "10"
        ]

        try:
            run_with_quiet_logging(cmd, operation_name="CPU search")

            # Stage 4: Convert results (CPU, ~5-30s)
            logging.info("📋 Stage 4/4: Converting results to M8 format (~5-30s)")
            logging.info("   Input: Alignment DB | Output: BLAST-like TSV | Hardware: CPU")
            cmd = [self.binary_path, "convertalis", query_db, self.base_db, result_db, result_m8]
            run_with_logging(cmd)

            return self._process_search_results(result_m8, target_sequence, query_db, tmp_dir)
        except subprocess.CalledProcessError:
            logging.warning("CPU search failed, returning minimal MSA")
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

    def _process_search_results(self, result_m8: str, target_sequence: str, query_db_path: str, tmp_dir: str) -> msa_tool.MsaToolResult:
        """Process search results to A3M."""
        import time
        stage_start = time.time()
        logging.info("🧬 Stage 5/5: Generating Multiple Sequence Alignment (~10-60s)")
        logging.info("   Input: M8 results + databases | Output: A3M MSA | Hardware: CPU")
        logging.info("   Algorithm: MSA construction with redundancy filtering")

        result_a3m_path = os.path.join(tmp_dir, "result.a3m")

        cmd = [
            self.binary_path, "result2msa", query_db_path, self.base_db, result_m8.replace('.m8', ''), result_a3m_path
        ]

        try:
            run_with_quiet_logging(cmd, operation_name="Generate MSA")
            with open(result_a3m_path) as f:
                a3m_content = f.read()

            # Validate and fix A3M format
            a3m_content = self._validate_and_fix_a3m(a3m_content, target_sequence)

            # Final validation: try to parse with AlphaFold3 to ensure it's valid
            try:
                from alphafold3.data.parsers import parse_fasta
                test_sequences, test_descriptions = parse_fasta(a3m_content)
                sequence_count = len([line for line in a3m_content.split('\n') if line.startswith('>')])
                logging.info("✅ MSA generation completed in %.1fs | Output: A3M format with %d sequences", time.time() - stage_start, sequence_count)
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

        except subprocess.CalledProcessError:
            logging.warning("MSA generation failed, returning minimal MSA")
            return msa_tool.MsaToolResult(
                a3m=f">query\n{target_sequence}\n",
                target_sequence=target_sequence,
                e_value=self.e_value
            )

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

            # Show parsed sequences after parsing
            logging.warning(f"📊 Parsed sequences: {len(sequences)} total")
            for i, (desc, seq) in enumerate(zip(descriptions, sequences)):
                logging.warning(f"  {i}: {desc} -> len={len(seq)}, content={repr(seq[:50])}...")
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
        logging.info("🔬 Starting MSA search for sequence (%d aa)", seq_len)
        logging.info("📋 Pipeline: FASTA → createdb → search → convertalis → result2msa → A3M")
        logging.info("⏱️  Estimated time: GPU=5-60s, CPU=5-60min | Hardware: %s",
                    "GPU+CPU" if self.use_gpu else "CPU only")
        
        with tempfile.TemporaryDirectory() as tmp_dir:
            query_path = os.path.join(tmp_dir, "query.fasta")
            result_m8 = os.path.join(tmp_dir, "result.m8")
            
            with open(query_path, "w") as f:
                f.write(f">query\n{target_sequence}\n")
            
            if self.use_gpu:
                logging.info("🎯 Using GPU acceleration (search stage)")
                result = self._gpu_search(query_path, result_m8, tmp_dir, target_sequence)
            else:
                logging.info("💻 Using CPU-only mode (all stages)")
                result = self._cpu_search(query_path, result_m8, tmp_dir, target_sequence)

            # Log total time
            total_time = time.time() - start_time
            logging.info("⏱️  Total MSA search completed in %.1fs", total_time)
            return result

