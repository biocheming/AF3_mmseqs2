#!/usr/bin/env python3
"""Test script to verify MMseqs2 capabilities detection."""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), 'src'))

from alphafold3.data.tools.mmseqs2 import MMseqs2

def test_capabilities():
    """Test MMseqs2 capability detection."""
    # Test with a dummy database path (doesn't need to exist for capability test)
    test_db = "/tmp/nonexistent.fasta"

    try:
        mmseqs2 = MMseqs2(
            binary_path="/usr/local/bin/mmseqs",  # Adjust path as needed
            database_path=test_db,
            n_cpu=4,
            sensitivity=4.0,
            use_gpu=False
        )

        print(f"MMseqs2 supports direct A3M output: {mmseqs2.supports_direct_a3m}")
        print(f"MMseqs2 supports MSA format mode 6: {mmseqs2.supports_msa_format_mode_6}")

        if mmseqs2.supports_direct_a3m:
            print("✓ Will use fast convertalis --format-output a3m mode")
        elif mmseqs2.supports_msa_format_mode_6:
            print("✓ Will use optimized result2msa --msa-format-mode 6 mode")
        else:
            print("✓ Will use standard result2msa mode for compatibility")

        return True

    except Exception as e:
        print(f"Error testing capabilities: {e}")
        return False

if __name__ == "__main__":
    success = test_capabilities()
    sys.exit(0 if success else 1)
