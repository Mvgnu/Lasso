#!/usr/bin/env python3
"""
Scan GBK genomes for alt-frame ORFs overlapping housekeeping gene loci,
score with the full ruleset, and summarize conservation across genomes.
License: MIT
Author: Magnus Ohle
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import logging
import re
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from lasso_workbench.core.prediction import CorePredictor, RuleEngine, load_active_ruleset
from lasso_workbench.pipeline.orf_extraction import chunk_orfs

logger = logging.getLogger(__name__)

DEFAULT_HOUSEKEEPING = [
    "rpoB",
    "rpoC",
    "gyrA",
    "gyrB",
    "recA",
    "tufA",
    "groEL",
    "dnaK",
    "rpsL",
    "rplB",
    "rplC",
]

DEFAULT_HOUSEKEEPING_SYNONYMS = {
    "rpoB": [
        "RNA polymerase beta subunit",
        "DNA-directed RNA polymerase subunit beta",
    ],
    "rpoC": [
        "RNA polymerase beta' subunit",
        "DNA-directed RNA polymerase subunit beta'",
    ],
    "gyrA": ["DNA gyrase subunit A"],
    "gyrB": ["DNA gyrase subunit B"],
    "recA": ["recombinase A"],
    "tufA": ["elongation factor Tu", "EF-Tu"],
    "groEL": ["chaperonin GroEL", "Hsp60", "heat shock protein 60", "Cpn60"],
    "dnaK": ["chaperone protein DnaK", "Hsp70"],
    "rpsL": ["ribosomal protein S12"],
    "rplB": ["ribosomal protein L2"],
    "rplC": ["ribosomal protein L3"],
}
DEFAULT_GEOM_BINS = 20
MAX_EXAMPLE_GENOMES = 5


@dataclass(frozen=True)
class HousekeepingGene:
    name: str
    start: int
    end: int
    strand: str
    frame_id: int
    matched_field: str
    matched_text: str

def _normalize_gene_list(names: Iterable[str]) -> List[str]:
    return [name.strip() for name in names if name and name.strip()]


def _load_housekeeping_names(path: Optional[Path], extra: Sequence[str]) -> List[str]:
    names: List[str] = []
    if path:
        with path.open() as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                names.append(line)
    names.extend(extra)
    names = _normalize_gene_list(names)
    return names or DEFAULT_HOUSEKEEPING


def _expand_with_synonyms(names: Sequence[str], use_synonyms: bool) -> List[Tuple[str, str]]:
    expanded: List[Tuple[str, str]] = []
    for name in names:
        expanded.append((name, name))
        if use_synonyms:
            for alias in DEFAULT_HOUSEKEEPING_SYNONYMS.get(name, []):
                expanded.append((name, alias))
    return expanded


def _build_gene_matchers(aliases: Sequence[Tuple[str, str]], use_regex: bool) -> List[Tuple[str, re.Pattern]]:
    matchers: List[Tuple[str, re.Pattern]] = []
    for canonical, alias in aliases:
        pattern = alias if use_regex else re.escape(alias)
        matchers.append((canonical, re.compile(pattern, re.IGNORECASE)))
    return matchers


def _match_gene_name(
    fields: Sequence[Tuple[str, str]],
    matchers: Sequence[Tuple[str, re.Pattern]],
) -> Optional[Tuple[str, str, str]]:
    for gene_name, pattern in matchers:
        for field_name, field_value in fields:
            if field_value and pattern.search(field_value):
                return gene_name, field_name, field_value
    return None


def _collect_fields(qualifiers) -> dict[str, List[str]]:
    fields: dict[str, List[str]] = {}
    for key in ("gene", "product", "locus_tag", "protein_id", "note", "function", "gene_synonym", "db_xref"):
        for value in qualifiers.get(key, []):
            if value is None:
                continue
            fields.setdefault(key, []).append(str(value))
    return fields


def _select_gene_name(fields: dict[str, List[str]], mode: str) -> Optional[Tuple[str, str]]:
    if mode == "any":
        for key in ("gene", "locus_tag", "protein_id", "product"):
            values = fields.get(key)
            if values:
                return key, values[0].strip()
        return None
    values = fields.get(mode)
    if not values:
        return None
    return mode, values[0].strip()


def _gene_frame_id(start: int, end: int, strand: str, seq_len: int) -> Optional[int]:
    if strand == "+":
        return start % 3
    if strand == "-":
        frame = (seq_len - end) % 3
        return -(frame + 1)
    return None


def _extract_housekeeping_genes(
    record,
    matchers: Sequence[Tuple[str, re.Pattern]],
    use_gbk_names: bool,
    gene_field: str,
) -> List[HousekeepingGene]:
    genes: List[HousekeepingGene] = []
    seq_len = len(record.seq)
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
        if use_gbk_names:
            selected = _select_gene_name(fields, gene_field)
            if not selected:
                continue
            matched_field, matched_text = selected
            matched_name = matched_text
        else:
            flat_fields: List[Tuple[str, str]] = []
            for key, values in fields.items():
                for value in values:
                    flat_fields.append((key, value))
            matched = _match_gene_name(flat_fields, matchers)
            if not matched:
                continue
            matched_name, matched_field, matched_text = matched
        start = int(feat.location.start)
        end = int(feat.location.end)
        frame_id = _gene_frame_id(start, end, strand, seq_len)
        if frame_id is None:
            continue
        genes.append(
            HousekeepingGene(
                name=matched_name,
                start=start,
                end=end,
                strand=strand,
                frame_id=frame_id,
                matched_field=matched_field,
                matched_text=matched_text,
            )
        )
    genes.sort(key=lambda g: g.start)
    return genes


def _find_overlaps(
    genes: Sequence[HousekeepingGene],
    starts: Sequence[int],
    orf_start: int,
    orf_end: int,
) -> List[HousekeepingGene]:
    if not genes:
        return []
    # genes sorted by start; find candidates with start <= orf_end
    import bisect

    idx = bisect.bisect_right(starts, orf_end) - 1
    hits: List[HousekeepingGene] = []
    while idx >= 0:
        gene = genes[idx]
        if gene.end < orf_start:
            break
        if orf_start < gene.end and orf_end > gene.start:
            hits.append(gene)
        idx -= 1
    return hits


def _merge_windows(windows: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    if not windows:
        return []
    sorted_windows = sorted(windows, key=lambda w: w[0])
    merged = [sorted_windows[0]]
    for start, end in sorted_windows[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _build_gene_windows(
    genes: Sequence[HousekeepingGene],
    seq_len: int,
    max_aa: int,
    pad_nt: int,
) -> List[Tuple[int, int]]:
    expand = pad_nt if pad_nt > 0 else max_aa * 3
    windows = []
    for gene in genes:
        start = max(0, gene.start - expand)
        end = min(seq_len, gene.end + expand)
        windows.append((start, end))
    return _merge_windows(windows)


def _iter_orfs_in_windows(
    seq: Seq,
    windows: Sequence[Tuple[int, int]],
    min_aa: int,
    max_aa: int,
) -> Iterator:
    for window_start, window_end in windows:
        window_seq = seq[window_start:window_end]
        yield from chunk_orfs(window_seq, "+", window_start, window_end, min_aa, max_aa)
        rc_seq = window_seq.reverse_complement()
        yield from chunk_orfs(rc_seq, "-", window_start, window_end, min_aa, max_aa)


def _score_full_rules(
    seq: str,
    predictor: CorePredictor,
) -> Tuple[float, Optional[str], List[str]]:
    pred = predictor.predict(seq, top_n=1, allow_inverted=True)
    best = pred.best_prediction
    if best is None:
        return 0.0, None, []
    return float(best.score), best.orientation, list(best.reasons)


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:12]


def _write_metadata(path: Path, payload: dict) -> None:
    with path.open("w") as handle:
        json.dump(payload, handle, indent=2, sort_keys=True)


def _note_unique_count(
    counts: dict,
    last_seen: dict,
    key,
    record_uid: int,
) -> bool:
    if last_seen.get(key) == record_uid:
        return False
    counts[key] = counts.get(key, 0) + 1
    last_seen[key] = record_uid
    return True


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Check alt-frame ORF conservation across housekeeping gene loci using the full ruleset."
    )
    parser.add_argument(
        "--gbk-dir",
        type=Path,
        default=Path("data/antismash_lasso/gbk"),
        help="Folder containing GBK files",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/alt_frame_conservation"),
        help="Output directory for hits and summaries",
    )
    parser.add_argument("--min-aa", type=int, default=20, help="Minimum ORF length (aa)")
    parser.add_argument("--max-aa", type=int, default=120, help="Maximum ORF length (aa)")
    parser.add_argument("--min-score", type=float, default=8.0, help="Minimum rule score to keep")
    parser.add_argument("--housekeeping-file", type=Path, default=None, help="File with housekeeping gene names")
    parser.add_argument(
        "--gene",
        action="append",
        default=[],
        help="Housekeeping gene name (repeatable). Overrides defaults if provided.",
    )
    parser.add_argument(
        "--gene-regex",
        action="store_true",
        help="Treat housekeeping gene names as regex patterns",
    )
    parser.add_argument(
        "--use-gbk-genes",
        action="store_true",
        help="Use gene names straight from GBK CDS qualifiers (no list matching).",
    )
    parser.add_argument(
        "--gene-field",
        choices=["gene", "product", "locus_tag", "protein_id", "any"],
        default="gene",
        help="Which GBK qualifier to use when --use-gbk-genes is set",
    )
    parser.add_argument(
        "--prepass-only",
        action="store_true",
        help="Only scan GBKs for gene presence, skip ORF extraction and scoring",
    )
    parser.add_argument(
        "--min-gene-genomes",
        type=int,
        default=3,
        help="Minimum genomes required for allowlist when using --prepass-only",
    )
    parser.add_argument(
        "--allowed-gene-file",
        type=Path,
        default=None,
        help="Optional: only analyze genes listed here (one per line)",
    )
    parser.add_argument(
        "--allowed-gene-out",
        type=Path,
        default=None,
        help="Write allowlist of genes present in >= min-gene-genomes",
    )
    parser.add_argument(
        "--geom-bins",
        type=int,
        default=DEFAULT_GEOM_BINS,
        help="Number of bins for geometry summary (default: 20)",
    )
    parser.add_argument(
        "--no-synonyms",
        action="store_true",
        help="Disable built-in housekeeping gene synonyms",
    )
    parser.add_argument(
        "--gene-pad-nt",
        type=int,
        default=0,
        help="Extra nucleotide padding around gene windows (if set, replaces max_aa*3 expansion)",
    )
    parser.add_argument(
        "--scan-full-genome",
        action="store_true",
        help="Scan entire genome instead of windows around housekeeping genes",
    )
    parser.add_argument(
        "--write-gene-matches",
        action="store_true",
        help="Write matched housekeeping gene features to gene_matches.tsv",
    )
    parser.add_argument("--limit", type=int, default=None, help="Limit number of GBK files to process")
    parser.add_argument("--max-records", type=int, default=None, help="Limit records per GBK file")
    parser.add_argument(
        "--write-sequences",
        action="store_true",
        help="Include protein sequences in summary outputs (larger files)",
    )
    parser.add_argument("--log-every", type=int, default=25, help="Log progress every N files")
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
    )

    if not args.gbk_dir.exists():
        logger.error("GBK directory not found: %s", args.gbk_dir)
        return 2

    names = _load_housekeeping_names(args.housekeeping_file, args.gene)
    if args.use_gbk_genes:
        aliases = []
        matchers = []
        logger.info("Using GBK gene names (field=%s)", args.gene_field)
    else:
        aliases = _expand_with_synonyms(names, use_synonyms=not args.no_synonyms)
        matchers = _build_gene_matchers(aliases, args.gene_regex)
        logger.info("Housekeeping genes: %s", ", ".join(names))

    ruleset = load_active_ruleset()
    predictor = CorePredictor(rule_engine=RuleEngine(ruleset))

    args.output_dir.mkdir(parents=True, exist_ok=True)
    hits_path = args.output_dir / "alt_frame_hits.tsv"
    summary_genes_path = args.output_dir / "summary_genes.tsv"
    summary_gene_matches_path = args.output_dir / "summary_gene_matches.tsv"
    summary_sequences_path = args.output_dir / "summary_sequences.tsv"
    summary_genes_geom_path = args.output_dir / "summary_genes_geom.tsv"
    summary_sequences_geom_path = args.output_dir / "summary_sequences_geom.tsv"
    metadata_path = args.output_dir / "run_metadata.json"
    gene_matches_path = args.output_dir / "gene_matches.tsv"
    allowlist_path = args.allowed_gene_out or (
        args.output_dir / f"allowed_genes_min{args.min_gene_genomes}.txt"
        if args.prepass_only
        else None
    )

    gbk_files = sorted(args.gbk_dir.glob("*.gbk"))
    if args.limit is not None:
        gbk_files = gbk_files[: args.limit]

    total_records = 0
    record_uid = 0
    genomes_seen: set[str] = set()

    seq_unique_counts: dict[Tuple[str, str], int] = {}
    seq_last_seen: dict[Tuple[str, str], int] = {}
    seq_examples: dict[Tuple[str, str], List[str]] = {}
    seq_hits: dict[Tuple[str, str], int] = {}
    seq_len: dict[Tuple[str, str], int] = {}
    seq_text: dict[Tuple[str, str], str] = {}
    gene_hits: dict[str, int] = {}
    gene_unique_counts: dict[str, int] = {}
    gene_last_seen: dict[str, int] = {}
    geom_unique_counts: dict[Tuple, int] = {}
    geom_last_seen: dict[Tuple, int] = {}
    geom_examples: dict[Tuple, List[str]] = {}
    geom_hits: dict[Tuple, int] = {}
    geom_gene_hits: dict[str, int] = {}
    geom_gene_unique_counts: dict[str, int] = {}
    geom_gene_last_seen: dict[str, int] = {}
    gene_feature_counts: dict[str, int] = {}
    gene_feature_unique_counts: dict[str, int] = {}
    gene_feature_last_seen: dict[str, int] = {}

    start_time = time.time()

    allowed_genes: set[str] = set()
    if args.allowed_gene_file:
        if not args.allowed_gene_file.exists():
            logger.error("Allowed gene file not found: %s", args.allowed_gene_file)
            return 2
        with args.allowed_gene_file.open() as handle:
            for line in handle:
                line = line.strip()
                if line:
                    allowed_genes.add(line)
        logger.info("Loaded %d allowed genes", len(allowed_genes))

    gene_writer = None
    gene_handle = None
    if args.write_gene_matches:
        gene_handle = gene_matches_path.open("w", newline="")
        gene_writer = csv.writer(gene_handle, delimiter="\t")
        gene_writer.writerow(
            [
                "record_id",
                "gbk_file",
                "record_index",
                "gene_name",
                "gene_start",
                "gene_end",
                "gene_strand",
                "gene_frame",
                "matched_field",
                "matched_text",
            ]
        )

    with hits_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        header = [
            "record_id",
            "gbk_file",
            "record_index",
            "gene_name",
            "gene_start",
            "gene_end",
            "gene_strand",
            "gene_frame",
            "orf_start",
            "orf_end",
            "orf_strand",
            "orf_frame",
            "orf_aa_len",
            "orf_nt_len",
            "match_type",
            "rule_score",
            "rule_orientation",
            "seq_hash",
            "rule_reasons",
        ]
        if args.write_sequences:
            header.append("protein_sequence")
        writer.writerow(header)

        for file_idx, gbk_path in enumerate(gbk_files, start=1):
            if file_idx % max(1, args.log_every) == 0:
                logger.info("Processing %s (%d/%d)", gbk_path.name, file_idx, len(gbk_files))

            try:
                records = SeqIO.parse(str(gbk_path), "genbank")
            except Exception as exc:
                logger.warning("Skipping %s (parse error: %s)", gbk_path, exc)
                continue

            record_index = 0
            while True:
                try:
                    record = next(records)
                except StopIteration:
                    break
                except Exception as exc:
                    logger.warning(
                        "Skipping remainder of %s (parse error: %s)", gbk_path.name, exc
                    )
                    break
                record_index += 1
                if args.max_records and record_index > args.max_records:
                    break
                total_records += 1
                record_uid += 1
                record_id = record.id or f"{gbk_path.stem}_record{record_index}"
                genomes_seen.add(record_id)

                genes = _extract_housekeeping_genes(
                    record,
                    matchers,
                    use_gbk_names=args.use_gbk_genes,
                    gene_field=args.gene_field,
                )
                if not genes:
                    continue
                if args.prepass_only:
                    for gene in genes:
                        gene_feature_counts[gene.name] = gene_feature_counts.get(gene.name, 0) + 1
                        _note_unique_count(
                            gene_feature_unique_counts,
                            gene_feature_last_seen,
                            gene.name,
                            record_uid,
                        )
                        if gene_writer:
                            gene_writer.writerow(
                                [
                                    record_id,
                                    gbk_path.name,
                                    record_index,
                                    gene.name,
                                    gene.start,
                                    gene.end,
                                    gene.strand,
                                    gene.frame_id,
                                    gene.matched_field,
                                    gene.matched_text,
                                ]
                            )
                    continue

                if allowed_genes:
                    genes = [gene for gene in genes if gene.name in allowed_genes]
                    if not genes:
                        continue

                for gene in genes:
                    gene_feature_counts[gene.name] = gene_feature_counts.get(gene.name, 0) + 1
                    _note_unique_count(
                        gene_feature_unique_counts,
                        gene_feature_last_seen,
                        gene.name,
                        record_uid,
                    )
                    if gene_writer:
                        gene_writer.writerow(
                            [
                                record_id,
                                gbk_path.name,
                                record_index,
                                gene.name,
                                gene.start,
                                gene.end,
                                gene.strand,
                                gene.frame_id,
                                gene.matched_field,
                                gene.matched_text,
                            ]
                        )
                starts = [gene.start for gene in genes]

                seq = record.seq
                if args.scan_full_genome:
                    windows = [(0, len(seq))]
                else:
                    windows = _build_gene_windows(
                        genes,
                        len(seq),
                        max_aa=args.max_aa,
                        pad_nt=args.gene_pad_nt,
                    )
                for orf in _iter_orfs_in_windows(seq, windows, args.min_aa, args.max_aa):
                    overlaps = _find_overlaps(
                        genes,
                        starts,
                        int(orf.genomic_start),
                        int(orf.genomic_end),
                    )
                    if not overlaps:
                        continue

                    score, orientation, reasons = _score_full_rules(orf.protein, predictor)
                    if score < args.min_score:
                        continue

                    seq_hash = _seq_hash(orf.protein)
                    for gene in overlaps:
                        if int(orf.frame) == gene.frame_id:
                            continue
                        match_type = "antisense" if orf.strand != gene.strand else "out_of_frame"
                        key = (gene.name, seq_hash)
                        if _note_unique_count(seq_unique_counts, seq_last_seen, key, record_uid):
                            examples = seq_examples.setdefault(key, [])
                            if len(examples) < MAX_EXAMPLE_GENOMES:
                                examples.append(record_id)
                        seq_hits[key] = seq_hits.get(key, 0) + 1
                        seq_len.setdefault(key, len(orf.protein))
                        if args.write_sequences:
                            seq_text.setdefault(key, orf.protein)
                        gene_hits[gene.name] = gene_hits.get(gene.name, 0) + 1
                        _note_unique_count(gene_unique_counts, gene_last_seen, gene.name, record_uid)

                        gene_len = gene.end - gene.start
                        if gene_len > 0 and args.geom_bins > 0:
                            ov_start = max(gene.start, int(orf.genomic_start))
                            ov_end = min(gene.end, int(orf.genomic_end))
                            if ov_end > ov_start:
                                rel_start = (ov_start - gene.start) / gene_len
                                rel_end = (ov_end - gene.start) / gene_len
                                bin_start = int(rel_start * args.geom_bins)
                                bin_end = int(rel_end * args.geom_bins)
                                if bin_start < 0:
                                    bin_start = 0
                                if bin_end < 0:
                                    bin_end = 0
                                max_bin = args.geom_bins - 1
                                if bin_start > max_bin:
                                    bin_start = max_bin
                                if bin_end > max_bin:
                                    bin_end = max_bin
                                geom_key = (
                                    gene.name,
                                    seq_hash,
                                    match_type,
                                    gene.strand,
                                    gene.frame_id,
                                    orf.strand,
                                    int(orf.frame),
                                    bin_start,
                                    bin_end,
                                )
                                if _note_unique_count(
                                    geom_unique_counts, geom_last_seen, geom_key, record_uid
                                ):
                                    examples = geom_examples.setdefault(geom_key, [])
                                    if len(examples) < MAX_EXAMPLE_GENOMES:
                                        examples.append(record_id)
                                geom_hits[geom_key] = geom_hits.get(geom_key, 0) + 1
                                geom_gene_hits[gene.name] = geom_gene_hits.get(gene.name, 0) + 1
                                _note_unique_count(
                                    geom_gene_unique_counts,
                                    geom_gene_last_seen,
                                    gene.name,
                                    record_uid,
                                )

                        row = [
                            record_id,
                            gbk_path.name,
                            record_index,
                            gene.name,
                            gene.start,
                            gene.end,
                            gene.strand,
                            gene.frame_id,
                            int(orf.genomic_start),
                            int(orf.genomic_end),
                            orf.strand,
                            int(orf.frame),
                            int(orf.aa_len),
                            int(orf.nt_len),
                            match_type,
                            f"{score:.2f}",
                            orientation or "",
                            seq_hash,
                            ";".join(reasons),
                        ]
                        if args.write_sequences:
                            row.append(orf.protein)
                        writer.writerow(row)

    if gene_handle:
        gene_handle.close()

    summary_rows = []
    for gene_name in sorted(gene_unique_counts):
        genomes_with_hit = gene_unique_counts.get(gene_name, 0)
        total_hits = gene_hits.get(gene_name, 0)
        unique_sequences = sum(1 for key in seq_hits if key[0] == gene_name)
        conserved_sequences = sum(
            1
            for key, count in seq_unique_counts.items()
            if key[0] == gene_name and count >= 2
        )
        summary_rows.append(
            {
                "gene_name": gene_name,
                "genomes_with_hit": genomes_with_hit,
                "total_hits": total_hits,
                "unique_sequences": unique_sequences,
                "conserved_sequences": conserved_sequences,
            }
        )
    summary_columns = [
        "gene_name",
        "genomes_with_hit",
        "total_hits",
        "unique_sequences",
        "conserved_sequences",
    ]
    pd.DataFrame(summary_rows, columns=summary_columns).to_csv(summary_genes_path, sep="\t", index=False)

    gene_match_rows = []
    for gene_name in sorted(gene_feature_unique_counts):
        gene_match_rows.append(
            {
                "gene_name": gene_name,
                "records_with_gene": gene_feature_unique_counts.get(gene_name, 0),
                "gene_feature_count": gene_feature_counts.get(gene_name, 0),
            }
        )
    pd.DataFrame(
        gene_match_rows,
        columns=["gene_name", "records_with_gene", "gene_feature_count"],
    ).to_csv(summary_gene_matches_path, sep="\t", index=False)

    if args.prepass_only and allowlist_path:
        allowed = sorted(
            [
                gene
                for gene, count in gene_feature_unique_counts.items()
                if count >= args.min_gene_genomes
            ]
        )
        with allowlist_path.open("w") as handle:
            for gene in allowed:
                handle.write(f"{gene}\n")
        logger.info("Wrote %s", allowlist_path)

    seq_rows = []
    for (gene_name, seq_hash), count in seq_unique_counts.items():
        row = {
            "gene_name": gene_name,
            "seq_hash": seq_hash,
            "aa_length": seq_len.get((gene_name, seq_hash), 0),
            "genomes": count,
            "hits": seq_hits.get((gene_name, seq_hash), 0),
            "conserved": count >= 2,
            "example_genomes": ",".join(seq_examples.get((gene_name, seq_hash), [])),
        }
        if args.write_sequences:
            row["protein_sequence"] = seq_text.get((gene_name, seq_hash), "")
        seq_rows.append(row)
    seq_columns = [
        "gene_name",
        "seq_hash",
        "aa_length",
        "genomes",
        "hits",
        "conserved",
        "example_genomes",
    ]
    if args.write_sequences:
        seq_columns.append("protein_sequence")
    seq_df = pd.DataFrame(seq_rows, columns=seq_columns)
    if not seq_df.empty:
        seq_df = seq_df.sort_values(
            ["gene_name", "genomes", "hits"], ascending=[True, False, False]
        )
    seq_df.to_csv(summary_sequences_path, sep="\t", index=False)

    geom_rows = []
    for geom_key, count in geom_unique_counts.items():
        (
            gene_name,
            seq_hash,
            match_type,
            gene_strand,
            gene_frame,
            orf_strand,
            orf_frame,
            bin_start,
            bin_end,
        ) = geom_key
        geom_rows.append(
            {
                "gene_name": gene_name,
                "seq_hash": seq_hash,
                "match_type": match_type,
                "gene_strand": gene_strand,
                "gene_frame": gene_frame,
                "orf_strand": orf_strand,
                "orf_frame": orf_frame,
                "bin_start": bin_start,
                "bin_end": bin_end,
                "genomes": count,
                "hits": geom_hits.get(geom_key, 0),
                "conserved": count >= 2,
                "example_genomes": ",".join(geom_examples.get(geom_key, [])),
            }
        )
    geom_columns = [
        "gene_name",
        "seq_hash",
        "match_type",
        "gene_strand",
        "gene_frame",
        "orf_strand",
        "orf_frame",
        "bin_start",
        "bin_end",
        "genomes",
        "hits",
        "conserved",
        "example_genomes",
    ]
    geom_df = pd.DataFrame(geom_rows, columns=geom_columns)
    if not geom_df.empty:
        geom_df = geom_df.sort_values(
            ["gene_name", "genomes", "hits"], ascending=[True, False, False]
        )
    geom_df.to_csv(summary_sequences_geom_path, sep="\t", index=False)

    geom_gene_rows = []
    for gene_name in sorted(geom_gene_unique_counts):
        genomes_with_hit = geom_gene_unique_counts.get(gene_name, 0)
        total_hits = geom_gene_hits.get(gene_name, 0)
        unique_sequences = sum(1 for key in geom_hits if key[0] == gene_name)
        conserved_sequences = sum(
            1
            for key, count in geom_unique_counts.items()
            if key[0] == gene_name and count >= 2
        )
        geom_gene_rows.append(
            {
                "gene_name": gene_name,
                "genomes_with_hit": genomes_with_hit,
                "total_hits": total_hits,
                "unique_sequences": unique_sequences,
                "conserved_sequences": conserved_sequences,
            }
        )
    pd.DataFrame(
        geom_gene_rows,
        columns=[
            "gene_name",
            "genomes_with_hit",
            "total_hits",
            "unique_sequences",
            "conserved_sequences",
        ],
    ).to_csv(summary_genes_geom_path, sep="\t", index=False)

    elapsed = time.time() - start_time
    _write_metadata(
        metadata_path,
        {
            "gbk_dir": str(args.gbk_dir),
            "gbk_files": len(gbk_files),
            "total_records": total_records,
            "genomes_seen": len(genomes_seen),
            "min_aa": args.min_aa,
            "max_aa": args.max_aa,
            "min_score": args.min_score,
            "housekeeping_genes": names,
            "housekeeping_synonyms": not args.no_synonyms,
            "use_gbk_genes": args.use_gbk_genes,
            "gene_field": args.gene_field,
            "gene_pad_nt": args.gene_pad_nt,
            "scan_full_genome": args.scan_full_genome,
            "prepass_only": args.prepass_only,
            "min_gene_genomes": args.min_gene_genomes,
            "allowed_gene_file": str(args.allowed_gene_file) if args.allowed_gene_file else None,
            "geom_bins": args.geom_bins,
            "ruleset_version": ruleset.version,
            "ruleset_description": ruleset.description,
            "ruleset_enabled_rules": [
                name for name, rule in ruleset.rules.items() if rule.enabled
            ],
            "runtime_seconds": round(elapsed, 2),
        },
    )

    logger.info("Wrote %s", hits_path)
    logger.info("Wrote %s", summary_genes_path)
    logger.info("Wrote %s", summary_gene_matches_path)
    logger.info("Wrote %s", summary_sequences_path)
    logger.info("Wrote %s", summary_genes_geom_path)
    logger.info("Wrote %s", summary_sequences_geom_path)
    logger.info("Wrote %s", metadata_path)
    logger.info("Done in %.1fs", elapsed)
    return 0


if __name__ == "__main__":
    sys.exit(main())
