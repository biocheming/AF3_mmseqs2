#!/usr/bin/env python3
"""
Test script to verify fast A3M generation is working correctly.
"""

import os
import sys
import tempfile

# Add project path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from alphafold3.data.tools.mmseqs2 import MMseqs2

def test_fast_a3m():
    """Test fast A3M generation with real database."""

    # Test sequence
    test_sequence = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"

    # Use the test database
    db_path = "src/alphafold3/test_data/miniature_databases/uniref90__subsampled_1000.fasta"

    if not os.path.exists(db_path):
        print(f"❌ Test database not found: {db_path}")
        return False

    print("🧪 Testing fast A3M generation...")
    print(f"   Sequence length: {len(test_sequence)}")
    print(f"   Database: {db_path}")

    try:
        # Create MMseqs2 instance
        mmseqs2 = MMseqs2(
            binary_path="/home/ldl/mmseqs/bin/mmseqs",
            database_path=db_path,
            n_cpu=4,
            e_value=0.0001,
            max_sequences=10000,
            sensitivity=4.0,
            gpu_devices=None,
            use_gpu=False
        )

        # Run MSA generation
        result = mmseqs2.query(test_sequence)

        # Check result
        if result and result.a3m:
            sequence_count = len([line for line in result.a3m.split('\n') if line.startswith('>')])
            print(f"✅ Fast A3M generation successful!")
            print(f"   Sequences in MSA: {sequence_count}")
            print(f"   A3M size: {len(result.a3m)} characters")

            # Check if it contains the query sequence
            if test_sequence in result.a3m:
                print("✅ Query sequence found in A3M")
            else:
                print("⚠️  Query sequence not found in A3M")

            return True
        else:
            print("❌ No A3M result returned")
            return False

    except Exception as e:
        print(f"❌ Test failed with error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = test_fast_a3m()
    sys.exit(0 if success else 1)
