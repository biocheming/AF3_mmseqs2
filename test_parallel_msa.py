#!/usr/bin/env python3
"""Test script for parallel MSA generation across multiple databases (GPT Scheme 3)."""

import sys
import os
import time
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from alphafold3.data.tools.mmseqs2 import MMseqs2

def test_parallel_msa():
    """Test parallel MSA generation across multiple databases."""
    # Test sequence
    target_sequence = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    # Example database configurations
    # In practice, these would come from your actual database paths
    database_configs = [
        {
            'binary_path': '/usr/local/bin/mmseqs',
            'database_path': '/opt/data2/af3_db/uniref90_2022_05.fa',  # This would be actual path
            'n_cpu': 8,
            'e_value': 0.0001,
            'max_sequences': 3000,
            'sensitivity': 7.5,
            'gpu_devices': ['0'],  # GPU enabled
            'use_gpu': True
        },
        {
            'binary_path': '/usr/local/bin/mmseqs',
            'database_path': '/opt/data2/af3_db/bfd-first_non_consensus_sequences.fasta',
            'n_cpu': 6,
            'e_value': 0.0001,
            'max_sequences': 2000,
            'sensitivity': 7.5,
            'gpu_devices': None,  # CPU only
            'use_gpu': False
        },
        {
            'binary_path': '/usr/local/bin/mmseqs',
            'database_path': '/opt/data2/af3_db/mgy_clusters_2022_05.fa',
            'n_cpu': 4,
            'e_value': 0.0001,
            'max_sequences': 1500,
            'sensitivity': 7.5,
            'gpu_devices': ['1'],  # Different GPU
            'use_gpu': True
        }
    ]

    print("🚀 Testing parallel MSA generation across multiple databases...")
    print(f"Target sequence length: {len(target_sequence)}")
    print(f"Number of databases: {len(database_configs)}")

    start_time = time.time()

    try:
        # Use the new parallel method
        results = MMseqs2.query_multiple_databases_parallel(target_sequence, database_configs)

        total_time = time.time() - start_time

        print(f"\n✅ Parallel MSA generation completed in {total_time:.1f}s")
        print("\nResults summary:")

        total_sequences = 0
        for db_name, result in results.items():
            sequence_count = len([line for line in result.a3m.split('\n') if line.startswith('>')])
            total_sequences += sequence_count
            print(f"  {db_name}: {sequence_count} sequences")

        print(f"\nTotal sequences across all databases: {total_sequences}")
        print(f"Average time per database: {total_time / len(database_configs):.1f}s")

        return True

    except Exception as e:
        print(f"❌ Parallel MSA test failed: {e}")
        return False

if __name__ == "__main__":
    # Note: This test will fail because the database paths don't exist
    # It's just for demonstrating the API
    print("⚠️  This is a demonstration script. Actual database paths need to be configured.")
    print("    The script shows how to use the new parallel MSA generation feature.")

    success = test_parallel_msa()
    sys.exit(0 if success else 1)
