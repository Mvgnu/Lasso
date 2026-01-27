"""
Utility modules for lasso peptide analysis.

Modules:
- translation: Bacterial translation with proper start codon handling
"""

from lasso_workbench.utils.translation import (
    translate_bacterial,
    START_CODONS,
    STOP_CODONS,
    BACTERIAL_CODON_TABLE,
)
from lasso_workbench.utils.sequence_io import (
    read_fasta_pairs,
    write_fasta_pairs,
    write_fasta_entries,
    read_genbank_record,
    iter_genbank_records,
)

__all__ = [
    "translate_bacterial",
    "START_CODONS",
    "STOP_CODONS",
    "BACTERIAL_CODON_TABLE",
    "read_fasta_pairs",
    "write_fasta_pairs",
    "write_fasta_entries",
    "read_genbank_record",
    "iter_genbank_records",
]
