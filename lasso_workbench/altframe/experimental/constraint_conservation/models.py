from __future__ import annotations

from dataclasses import dataclass
from typing import List


@dataclass
class GeneInstance:
    record_uid: str
    record_id: str
    gbk_file: str
    gene_name: str
    gene_start: int
    gene_end: int
    gene_strand: str
    gene_len: int
    gene_genomic: str  # genomic orientation, length = gene_len
    cds_prefix: str
    codons: List[str]  # CDS orientation
    aas: List[str]
    remainder: str


@dataclass
class WindowInstance:
    rel_start: int
    rel_end: int
    gene_strand: str
    cds_prefix: str
    codons: List[str]
    aas: List[str]
    remainder: str
    cds_len: int
    window_seq: str


@dataclass(frozen=True)
class LocusKey:
    gene_name: str
    match_type: str
    orf_strand: str
    orf_frame: int
    bin_start: int
    bin_end: int
