#!/usr/bin/env python3
"""
Constraint-aware conservation test for alt-frame candidates (locus-level).

Defines candidates by locus (gene + bin range + alt frame), not by seq_hash.
For each locus:
- Extract alt-frame peptides from all genomes with the gene (windowed by bins).
- Compute observed ORF survival (no stop) + mean pairwise identity among survivors.
- Generate nulls by synonymous codon permutation within each genome (preserve AA + codon counts per gene).
- Compare observed to null distributions (p-values, z-score).
"""
from __future__ import annotations

import argparse
import csv
import logging
import math
import random
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Data import CodonTable

from lasso_workbench.utils.translation import translate_bacterial

logger = logging.getLogger(__name__)


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
    codons: List[str]  # CDS orientation
    aas: List[str]
    remainder: str


@dataclass(frozen=True)
class LocusKey:
    gene_name: str
    match_type: str
    orf_strand: str
    orf_frame: int
    bin_start: int
    bin_end: int


def _reverse_complement(seq: str) -> str:
    comp = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
    return seq.translate(comp)[::-1]


def _mean_pairwise_identity(seqs: Sequence[str]) -> float:
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


def _collect_fields(qualifiers) -> dict[str, List[str]]:
    fields: dict[str, List[str]] = {}
    for key in ("gene", "product", "locus_tag", "protein_id", "note", "function", "gene_synonym", "db_xref"):
        for value in qualifiers.get(key, []):
            if value is None:
                continue
            fields.setdefault(key, []).append(str(value))
    return fields


def _select_gene_name(fields: dict[str, List[str]], mode: str) -> Optional[str]:
    if mode == "any":
        for key in ("gene", "locus_tag", "protein_id", "product"):
            values = fields.get(key)
            if values:
                return values[0].strip()
        return None
    values = fields.get(mode)
    if not values:
        return None
    return values[0].strip()


def _compute_bins(
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
    rel_start = (ov_start - gene_start) / gene_len
    rel_end = (ov_end - gene_start) / gene_len
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


def _window_from_bins(gene_start: int, gene_end: int, bin_start: int, bin_end: int, bins: int) -> Tuple[int, int]:
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


def _translate_window(window_seq: str, orf_strand: str, orf_frame: int) -> str:
    seq = window_seq
    if orf_strand == "-":
        seq = _reverse_complement(seq)
        frame_offset = abs(int(orf_frame)) - 1
    else:
        if orf_frame < 0:
            frame_offset = abs(int(orf_frame)) - 1
        else:
            frame_offset = int(orf_frame)
    frame_offset = frame_offset % 3
    if frame_offset:
        seq = seq[frame_offset:]
    if len(seq) < 3:
        return ""
    return translate_bacterial(seq, start_codon_to_met=True)


def _synonym_maps() -> Tuple[Dict[str, List[str]], Dict[str, str]]:
    table = CodonTable.unambiguous_dna_by_id[11]
    codon_to_aa = {codon.upper(): aa for codon, aa in table.forward_table.items()}
    aa_to_codons: Dict[str, List[str]] = defaultdict(list)
    for codon, aa in codon_to_aa.items():
        aa_to_codons[aa].append(codon)
    return aa_to_codons, codon_to_aa


def _shuffle_synonymous_codons(codons: Sequence[str], aas: Sequence[str], rng: random.Random) -> List[str]:
    aa_to_positions: Dict[str, List[int]] = defaultdict(list)
    aa_to_codons: Dict[str, List[str]] = defaultdict(list)
    for idx, (codon, aa) in enumerate(zip(codons, aas)):
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


def _p_value(count_ge: int, n: int) -> float:
    if n <= 0:
        return float("nan")
    return (count_ge + 1) / (n + 1)


def _scan_hits_for_loci(
    hits_path: Path,
    geom_bins: int,
) -> Tuple[Dict[LocusKey, int], Dict[LocusKey, int]]:
    unique_counts: Dict[LocusKey, int] = {}
    last_seen: Dict[LocusKey, str] = {}
    hit_counts: Dict[LocusKey, int] = defaultdict(int)

    with hits_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            try:
                gene_start = int(row["gene_start"])
                gene_end = int(row["gene_end"])
                orf_start = int(row["orf_start"])
                orf_end = int(row["orf_end"])
            except (KeyError, ValueError):
                continue
            bins = _compute_bins(gene_start, gene_end, orf_start, orf_end, geom_bins)
            if not bins:
                continue
            bin_start, bin_end = bins
            key = LocusKey(
                gene_name=row["gene_name"],
                match_type=row["match_type"],
                orf_strand=row["orf_strand"],
                orf_frame=int(row["orf_frame"]),
                bin_start=bin_start,
                bin_end=bin_end,
            )
            record_uid = f"{row['gbk_file']}::{row['record_id']}"
            hit_counts[key] += 1
            if last_seen.get(key) != record_uid:
                last_seen[key] = record_uid
                unique_counts[key] = unique_counts.get(key, 0) + 1

    return unique_counts, hit_counts


def _select_top_loci(
    unique_counts: Dict[LocusKey, int],
    hit_counts: Dict[LocusKey, int],
    min_genomes: int,
    max_candidates: int,
) -> List[LocusKey]:
    rows: List[Tuple[int, int, LocusKey]] = []
    for key, genomes in unique_counts.items():
        if genomes < min_genomes:
            continue
        rows.append((genomes, hit_counts.get(key, 0), key))
    rows.sort(key=lambda x: (x[0], x[1]), reverse=True)
    if max_candidates and len(rows) > max_candidates:
        rows = rows[:max_candidates]
    return [key for _, _, key in rows]


def _extract_gene_instances(
    gbk_dir: Path,
    gene_names: set[str],
    gene_field: str,
) -> Dict[str, List[GeneInstance]]:
    _, codon_to_aa = _synonym_maps()
    instances: Dict[str, List[GeneInstance]] = defaultdict(list)
    seen: Dict[str, set[str]] = defaultdict(set)

    for gbk_path in sorted(gbk_dir.glob("*.gbk")):
        try:
            for record_index, record in enumerate(SeqIO.parse(str(gbk_path), "genbank"), start=1):
                record_id = record.id or f"{gbk_path.stem}_record{record_index}"
                record_uid = f"{gbk_path.name}::{record_id}"
                seq = str(record.seq).upper()
                for feat in record.features:
                    if feat.type != "CDS":
                        continue
                    strand_val = feat.location.strand
                    if strand_val == 1:
                        strand = "+"
                    elif strand_val == -1:
                        strand = "-"
                    else:
                        continue
                    qualifiers = feat.qualifiers or {}
                    fields = _collect_fields(qualifiers)
                    gene_name = _select_gene_name(fields, gene_field)
                    if not gene_name or gene_name not in gene_names:
                        continue
                    if record_uid in seen[gene_name]:
                        continue
                    start = int(feat.location.start)
                    end = int(feat.location.end)
                    if end <= start:
                        continue
                    gene_genomic = seq[start:end]
                    if not gene_genomic:
                        continue
                    cds_seq = gene_genomic if strand == "+" else _reverse_complement(gene_genomic)
                    codon_count = len(cds_seq) // 3
                    codons = [cds_seq[i:i + 3] for i in range(0, codon_count * 3, 3)]
                    aas = [codon_to_aa.get(codon, "X") for codon in codons]
                    remainder = cds_seq[codon_count * 3:]

                    instances[gene_name].append(
                        GeneInstance(
                            record_uid=record_uid,
                            record_id=record_id,
                            gbk_file=gbk_path.name,
                            gene_name=gene_name,
                            gene_start=start,
                            gene_end=end,
                            gene_strand=strand,
                            gene_len=end - start,
                            gene_genomic=gene_genomic,
                            codons=codons,
                            aas=aas,
                            remainder=remainder,
                        )
                    )
                    seen[gene_name].add(record_uid)
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed parsing %s: %s", gbk_path, exc)

    return instances


def _observed_for_locus(
    locus: LocusKey,
    gene_instances: Sequence[GeneInstance],
    geom_bins: int,
) -> Tuple[float, float, int, int, List[str]]:
    peptides: List[str] = []
    total = 0
    survived = 0

    for inst in gene_instances:
        window_start, window_end = _window_from_bins(
            inst.gene_start,
            inst.gene_end,
            locus.bin_start,
            locus.bin_end,
            geom_bins,
        )
        if window_end <= window_start:
            continue
        total += 1
        rel_start = window_start - inst.gene_start
        rel_end = window_end - inst.gene_start
        if rel_start < 0 or rel_end > inst.gene_len:
            continue
        window_seq = inst.gene_genomic[rel_start:rel_end]
        peptide = _translate_window(window_seq, locus.orf_strand, locus.orf_frame)
        if not peptide or "*" in peptide:
            continue
        survived += 1
        peptides.append(peptide)

    survival = survived / total if total else float("nan")
    identity = _mean_pairwise_identity(peptides)
    return survival, identity, total, survived, peptides


def _null_for_locus(
    locus: LocusKey,
    gene_instances: Sequence[GeneInstance],
    geom_bins: int,
    rng: random.Random,
    iterations: int,
) -> Tuple[List[float], List[float]]:
    survival_scores: List[float] = []
    identity_scores: List[float] = []

    for _ in range(iterations):
        peptides: List[str] = []
        total = 0
        survived = 0
        for inst in gene_instances:
            window_start, window_end = _window_from_bins(
                inst.gene_start,
                inst.gene_end,
                locus.bin_start,
                locus.bin_end,
                geom_bins,
            )
            if window_end <= window_start:
                continue
            shuffled_codons = _shuffle_synonymous_codons(inst.codons, inst.aas, rng)
            shuffled_cds = "".join(shuffled_codons) + inst.remainder
            if inst.gene_strand == "+":
                shuffled_genomic = shuffled_cds
            else:
                shuffled_genomic = _reverse_complement(shuffled_cds)
            total += 1
            rel_start = window_start - inst.gene_start
            rel_end = window_end - inst.gene_start
            if rel_start < 0 or rel_end > len(shuffled_genomic):
                continue
            window_seq = shuffled_genomic[rel_start:rel_end]
            peptide = _translate_window(window_seq, locus.orf_strand, locus.orf_frame)
            if not peptide or "*" in peptide:
                continue
            survived += 1
            peptides.append(peptide)

        survival = survived / total if total else float("nan")
        if not math.isnan(survival):
            survival_scores.append(survival)
        identity = _mean_pairwise_identity(peptides)
        if not math.isnan(identity):
            identity_scores.append(identity)

    return survival_scores, identity_scores


def main() -> int:
    parser = argparse.ArgumentParser(description="Constraint-aware conservation test for alt-frame loci.")
    parser.add_argument("--altframe-dir", type=Path, default=Path("results/alt_frame_conservation"))
    parser.add_argument("--min-genomes", type=int, default=20)
    parser.add_argument("--max-candidates", type=int, default=200)
    parser.add_argument("--geom-bins", type=int, default=20)
    parser.add_argument("--null-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--gbk-dir", type=Path, default=Path("data/antismash_lasso/gbk"))
    parser.add_argument("--hits-file", type=Path, default=None)
    parser.add_argument("--gene-field", choices=["gene", "locus_tag", "product", "any"], default="gene")
    parser.add_argument("--output-dir", type=Path, default=Path("results/altframe_constraint_tests"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    hits_path = args.hits_file or (args.altframe_dir / "alt_frame_hits.tsv")
    if not hits_path.exists():
        logger.error("Hits file not found: %s", hits_path)
        return 2

    unique_counts, hit_counts = _scan_hits_for_loci(hits_path, args.geom_bins)
    if not unique_counts:
        logger.error("No loci found in hits file.")
        return 2

    loci = _select_top_loci(unique_counts, hit_counts, args.min_genomes, args.max_candidates)
    if not loci:
        logger.error("No loci meet min-genomes filter.")
        return 2

    gene_names = {key.gene_name for key in loci}
    logger.info("Selected loci: %d | genes: %d", len(loci), len(gene_names))

    gene_instances = _extract_gene_instances(args.gbk_dir, gene_names, args.gene_field)
    if not gene_instances:
        logger.error("No gene instances found for selected loci.")
        return 2

    rng = random.Random(args.seed)

    rows: List[dict] = []
    for locus in loci:
        instances = gene_instances.get(locus.gene_name, [])
        if not instances:
            continue
        obs_survival, obs_identity, total, survived, _ = _observed_for_locus(
            locus,
            instances,
            args.geom_bins,
        )
        null_survival, null_identity = _null_for_locus(
            locus,
            instances,
            args.geom_bins,
            rng,
            args.null_iterations,
        )

        null_survival_mean = float(sum(null_survival) / len(null_survival)) if null_survival else float("nan")
        null_survival_std = float(
            math.sqrt(sum((x - null_survival_mean) ** 2 for x in null_survival) / len(null_survival))
        ) if null_survival else float("nan")
        null_identity_mean = float(sum(null_identity) / len(null_identity)) if null_identity else float("nan")
        null_identity_std = float(
            math.sqrt(sum((x - null_identity_mean) ** 2 for x in null_identity) / len(null_identity))
        ) if null_identity else float("nan")

        if null_identity_std and not math.isnan(obs_identity) and null_identity_std > 0:
            z_score = (obs_identity - null_identity_mean) / null_identity_std
        else:
            z_score = float("nan")

        if math.isnan(obs_survival):
            p_survival = float("nan")
        else:
            p_survival = _p_value(sum(1 for x in null_survival if x >= obs_survival), len(null_survival))
        if math.isnan(obs_identity):
            p_identity = float("nan")
        else:
            p_identity = _p_value(sum(1 for x in null_identity if x >= obs_identity), len(null_identity))

        rows.append({
            "gene_name": locus.gene_name,
            "match_type": locus.match_type,
            "orf_strand": locus.orf_strand,
            "orf_frame": locus.orf_frame,
            "bin_start": locus.bin_start,
            "bin_end": locus.bin_end,
            "genomes": unique_counts.get(locus, 0),
            "hits": hit_counts.get(locus, 0),
            "observed_n": total,
            "observed_survived": survived,
            "obs_survival": obs_survival,
            "obs_identity": obs_identity,
            "null_survival_mean": null_survival_mean,
            "null_survival_std": null_survival_std,
            "null_identity_mean": null_identity_mean,
            "null_identity_std": null_identity_std,
            "z_score": z_score,
            "p_survival": p_survival,
            "p_identity": p_identity,
            "null_samples_survival": len(null_survival),
            "null_samples_identity": len(null_identity),
        })

    if not rows:
        logger.error("No loci scored.")
        return 2

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "altframe_constraint_conservation.tsv"
    pd.DataFrame(rows).sort_values(by=["z_score"], ascending=False).to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote %d rows to %s", len(rows), out_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
