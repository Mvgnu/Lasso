from __future__ import annotations

import math
from typing import List, Sequence, Tuple

from Bio.Seq import Seq

from lasso_workbench.altframe.models import LocusKey, WindowInstance
from lasso_workbench.utils.translation import translate_bacterial


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement())


def mean_pairwise_identity(seqs: Sequence[str]) -> float:
    n = len(seqs)
    if n < 2:
        return float("nan")
    total = 0.0
    count = 0
    for i in range(n):
        s1 = seqs[i]
        if not s1:
            continue
        for j in range(i + 1, n):
            s2 = seqs[j]
            if not s2:
                continue
            m = min(len(s1), len(s2))
            if m == 0:
                continue
            matches = sum(1 for a, b in zip(s1[:m], s2[:m]) if a == b)
            total += matches / m
            count += 1
    if count == 0:
        return float("nan")
    return total / count


def _frame_offset(orf_frame: int, orf_strand: str, frame_mode: str) -> int:
    if orf_strand == "-":
        return (abs(int(orf_frame)) - 1) % 3
    if frame_mode == "one":
        return (abs(int(orf_frame)) - 1) % 3
    return int(orf_frame) % 3


def translate_window(window_seq: str, orf_strand: str, orf_frame: int, frame_mode: str) -> str:
    seq = window_seq
    if orf_strand == "-":
        seq = reverse_complement(seq)
    frame_offset = _frame_offset(orf_frame, orf_strand, frame_mode)
    if frame_offset:
        seq = seq[frame_offset:]
    if len(seq) < 3:
        return ""
    return translate_bacterial(seq, start_codon_to_met=True)


def observed_for_locus(
    locus: LocusKey,
    windows: Sequence[WindowInstance],
    frame_mode: str,
) -> Tuple[float, float, int, int, List[str]]:
    peptides: List[str] = []
    total = 0
    survived = 0

    for window in windows:
        total += 1
        peptide = translate_window(window.window_seq, locus.orf_strand, locus.orf_frame, frame_mode)
        if not peptide or "*" in peptide:
            continue
        survived += 1
        peptides.append(peptide)

    survival = survived / total if total else float("nan")
    identity = mean_pairwise_identity(peptides)
    return survival, identity, total, survived, peptides


def p_value(count_ge: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return (count_ge + 1) / (n + 1)


def summarize_null(scores: Sequence[float]) -> Tuple[float, float]:
    if not scores:
        return float("nan"), float("nan")
    mean = float(sum(scores) / len(scores))
    if len(scores) == 1:
        return mean, 0.0
    variance = sum((x - mean) ** 2 for x in scores) / len(scores)
    return mean, float(math.sqrt(variance))
