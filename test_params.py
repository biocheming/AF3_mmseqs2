#!/usr/bin/env python3
"""Test database-specific parameters."""

class MockMMseqs2:
    def __init__(self):
        self.max_sequences = 10000  # default
        self.sensitivity = 5.0

    def _get_database_specific_max_seqs(self, db_name: str) -> int:
        """Get database-specific max_sequences based on AlphaFold3 Jackhmmer settings."""
        db_name_lower = db_name.lower()

        # Database-specific max_sequences (balanced performance vs quality)
        if 'uniref90' in db_name_lower:
            return 5000   # Reduced from original 10,000 for speed
        elif 'mgy_clusters' in db_name_lower or 'mgnify' in db_name_lower:
            return 2000   # Keep conservative due to large database size
        elif 'bfd' in db_name_lower:
            return 3000   # Slightly increased from 2000
        elif 'uniprot' in db_name_lower:
            return 10000  # Reduced from original 50,000 for speed
        else:
            return self.max_sequences  # Fallback to default

# Test the logic
mmseqs = MockMMseqs2()
print('=== Database-specific max_sequences settings ===')
print(f'uniref90: {mmseqs._get_database_specific_max_seqs("uniref90_2022_05.fa")}')
print(f'mgy_clusters: {mmseqs._get_database_specific_max_seqs("mgy_clusters_2022_05.fa")}')
print(f'bfd-first: {mmseqs._get_database_specific_max_seqs("bfd-first_non_consensus_sequences.fasta")}')
print(f'uniprot: {mmseqs._get_database_specific_max_seqs("uniprot_all_2021_04.fa")}')
print(f'unknown: {mmseqs._get_database_specific_max_seqs("unknown_db.fa")}')
print(f'Default sensitivity: {mmseqs.sensitivity}')
print()
print('=== Comparison with original AlphaFold3 ===')
print('Original Jackhmmer settings:')
print('uniref90:     10,000 sequences')
print('mgnify:       5,000 sequences')
print('small_bfd:    5,000 sequences')
print('uniprot:      50,000 sequences')
print()
print('Our balanced settings:')
print('uniref90:     5,000 sequences (50% of original)')
print('mgy_clusters: 2,000 sequences (40% of original mgnify)')
print('bfd-first:    3,000 sequences (60% of original small_bfd)')
print('uniprot:      10,000 sequences (20% of original)')
print(f'Total: ~20,000 sequences (vs original ~70,000)')
print(f'Sensitivity: {mmseqs.sensitivity} (vs original Jackhmmer implicit high sensitivity)')
