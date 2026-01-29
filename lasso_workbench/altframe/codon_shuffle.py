from __future__ import annotations

import math
import random
from collections import defaultdict
from typing import Dict, List, Sequence, Tuple

from Bio.Data import CodonTable

from lasso_workbench.altframe.conservation import mean_pairwise_identity, reverse_complement, translate_window
from lasso_workbench.altframe.models import LocusKey, WindowInstance
from lasso_workbench.utils.translation import BACTERIAL_CODON_TABLE


def synonym_maps() -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    table = CodonTable.unambiguous_dna_by_id[BACTERIAL_CODON_TABLE]
    codon_to_aa = {codon.upper(): aa for codon, aa in table.forward_table.items()}
    aa_to_codons: Dict[str, List[str]] = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        aa_to_codons[aa].append(codon)
    return aa_to_codons, codon_to_aa


def shuffle_synonymous_codons(codons: Sequence[str], aas: Sequence[str], rng: random.Random) -> List[str]:
    aa_to_positions: Dict[str, List[int]] = defaultdict(list)
    aa_to_codons: Dict[str, List[str]] = defaultdict(list)
    for idx, (codon, aa) in enumerate(zip(codons, aas)):
        if idx == 0:
            continue
        if aa == "X":
            continue
        aa_to_positions[aa].append(idx)
        aa_to_codons[aa].append(codon)

    shuffled = list(codons)
    for aa, positions in aa_to_positions.items():
        codon_list = aa_to_codons[aa]
        shuffled_codons = codon_list[:]
        rng.shuffle(shuffled_codons)
        for pos, codon in zip(positions, shuffled_codons):
            shuffled[pos] = codon

    return shuffled


def null_for_locus(
    locus: LocusKey,
    windows: Sequence[WindowInstance],
    rng: random.Random,
    iterations: int,
    frame_mode: str,
) -> Tuple[List[float], List[float]]:
    survival_scores: List[float] = []
    identity_scores: List[float] = []

    for _ in range(iterations):
        peptides: List[str] = []
        total = 0
        survived = 0
        for window in windows:
            total += 1
            shuffled_codons = shuffle_synonymous_codons(window.codons, window.aas, rng)
            shuffled_cds = "".join(shuffled_codons) + window.remainder
            if window.gene_strand == "+":
                window_seq = shuffled_cds[window.rel_start:window.rel_end]
            else:
                seg = shuffled_cds[window.cds_len - window.rel_end : window.cds_len - window.rel_start]
                window_seq = reverse_complement(seg)
            peptide = translate_window(window_seq, locus.orf_strand, locus.orf_frame, frame_mode)
            if not peptide or "*" in peptide:
                continue
            survived += 1
            peptides.append(peptide)

        survival = survived / total if total else float("nan")
        if not math.isnan(survival):
            survival_scores.append(survival)

        identity = mean_pairwise_identity(peptides)
        if not math.isnan(identity):
            identity_scores.append(identity)

    return survival_scores, identity_scores
