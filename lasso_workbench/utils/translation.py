#!/usr/bin/env python3
"""Bacterial translation helpers backed by Biopython translation."""
from __future__ import annotations

from Bio.Seq import Seq

# Biopython translation table 11 = Bacterial/Archaeal/Plant Plastid genetic code.
BACTERIAL_CODON_TABLE = 11
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}
ALT_START_CODONS = {"GTG", "TTG"}


def is_alt_start_codon(codon: str) -> bool:
    return codon.upper() in ALT_START_CODONS


def translate_bacterial(
    dna: str | Seq,
    start_codon_to_met: bool = True,
    table: int = BACTERIAL_CODON_TABLE,
) -> str:
    """Translate DNA using bacterial codon table and normalize alt start codons to M."""
    dna_str = str(dna).upper()
    if len(dna_str) < 3:
        return ""

    protein = str(Seq(dna_str).translate(table=table))

    if protein.endswith("*"):
        protein = protein[:-1]

    if start_codon_to_met and protein and is_alt_start_codon(dna_str[:3]):
        protein = "M" + protein[1:]

    return protein
