from __future__ import annotations

from typing import List, Optional, Sequence, Tuple

from lasso_workbench.altframe.models import GeneInstance, LocusKey, WindowInstance


def compute_bins(
    gene_start: int,
    gene_end: int,
    orf_start: int,
    orf_end: int,
    bins: int,
) -> Optional[Tuple[int, int]]:
    gene_len = gene_end - gene_start
    if gene_len <= 0 or bins <= 0:
        return None
    ov_start = max(gene_start, orf_start)
    ov_end = min(gene_end, orf_end)
    if ov_end <= ov_start:
        return None
    # Coordinates are 0-based, end-exclusive. Use the last covered base for bin_end.
    last = ov_end - 1
    rel_start = (ov_start - gene_start) / gene_len
    rel_end = (last - gene_start) / gene_len
    bin_start = int(rel_start * bins)
    bin_end = int(rel_end * bins)
    if bin_start < 0:
        bin_start = 0
    if bin_end < 0:
        bin_end = 0
    max_bin = bins - 1
    if bin_start > max_bin:
        bin_start = max_bin
    if bin_end > max_bin:
        bin_end = max_bin
    return bin_start, bin_end


def window_from_bins(gene_start: int, gene_end: int, bin_start: int, bin_end: int, bins: int) -> Tuple[int, int]:
    gene_len = gene_end - gene_start
    if gene_len <= 0:
        return gene_start, gene_start
    start_frac = bin_start / bins
    end_frac = (bin_end + 1) / bins
    start = gene_start + int(start_frac * gene_len)
    end = gene_start + int(end_frac * gene_len)
    if end <= start:
        return gene_start, gene_start
    if start < gene_start:
        start = gene_start
    if end > gene_end:
        end = gene_end
    return start, end


def prepare_locus_windows(
    locus: LocusKey,
    gene_instances: Sequence[GeneInstance],
    geom_bins: int,
) -> List[WindowInstance]:
    windows: List[WindowInstance] = []
    for inst in gene_instances:
        window_start, window_end = window_from_bins(
            inst.gene_start,
            inst.gene_end,
            locus.bin_start,
            locus.bin_end,
            geom_bins,
        )
        if window_end <= window_start:
            continue
        rel_start = window_start - inst.gene_start
        rel_end = window_end - inst.gene_start
        if rel_start < 0 or rel_end > inst.gene_len:
            continue
        window_seq = inst.gene_genomic[rel_start:rel_end]
        windows.append(
            WindowInstance(
                rel_start=rel_start,
                rel_end=rel_end,
                gene_strand=inst.gene_strand,
                cds_prefix=inst.cds_prefix,
                codons=inst.codons,
                aas=inst.aas,
                remainder=inst.remainder,
                cds_len=inst.gene_len,
                window_seq=window_seq,
            )
        )
    return windows
