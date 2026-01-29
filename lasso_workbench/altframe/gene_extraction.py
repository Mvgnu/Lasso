from __future__ import annotations

import csv
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq

from lasso_workbench.altframe.binning import compute_bins
from lasso_workbench.altframe.codon_shuffle import synonym_maps
from lasso_workbench.altframe.models import GeneInstance, LocusKey

logger = logging.getLogger(__name__)


def collect_fields(qualifiers) -> dict[str, List[str]]:
    fields: dict[str, List[str]] = {}
    for key in ("gene", "product", "locus_tag", "protein_id", "note", "function", "gene_synonym", "db_xref"):
        for value in qualifiers.get(key, []):
            if value is None:
                continue
            fields.setdefault(key, []).append(str(value))
    return fields


def select_gene_name(fields: dict[str, List[str]], mode: str) -> Optional[str]:
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


def scan_hits_for_loci(
    hits_path: Path,
    geom_bins: int,
) -> Tuple[Dict[LocusKey, int], Dict[LocusKey, int], set[int]]:
    seen_uids: Dict[LocusKey, set[str]] = defaultdict(set)
    hit_counts: Dict[LocusKey, int] = defaultdict(int)
    frame_values: set[int] = set()

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
            bins = compute_bins(gene_start, gene_end, orf_start, orf_end, geom_bins)
            if not bins:
                continue
            bin_start, bin_end = bins
            orf_frame = int(row["orf_frame"])
            key = LocusKey(
                gene_name=row["gene_name"],
                match_type=row["match_type"],
                orf_strand=row["orf_strand"],
                orf_frame=orf_frame,
                bin_start=bin_start,
                bin_end=bin_end,
            )
            record_uid = f"{row['gbk_file']}::{row['record_id']}"
            hit_counts[key] += 1
            seen_uids[key].add(record_uid)
            frame_values.add(orf_frame)

    unique_counts = {k: len(v) for k, v in seen_uids.items()}
    return unique_counts, hit_counts, frame_values


def select_top_loci(
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


def extract_gene_instances(
    gbk_dir: Path,
    gene_names: set[str],
    gene_field: str,
) -> Dict[str, List[GeneInstance]]:
    _, codon_to_aa = synonym_maps()
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
                    fields = collect_fields(qualifiers)
                    gene_name = select_gene_name(fields, gene_field)
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
                    cds_seq = gene_genomic if strand == "+" else str(Seq(gene_genomic).reverse_complement())
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
