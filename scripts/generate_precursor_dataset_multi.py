#!/usr/bin/env python3
"""Generate multi-candidate precursor datasets anchored on lab core + peptidase loci.

Workflow:
1) Parse lab core table (semicolon CSV).
2) Locate peptidase loci in genomes by protein_id accession.
3) Rank loci per case (strict core matches + distance); keep all tied ambiguous loci.
4) Enumerate ORFs in ±window for retained loci, collect ORFs containing the core.
5) Emit strict validated FASTA (multi-candidate per core) + locus TSV + case summary TSV.

Matching rule:
  A (primary): I/L equivalence only, core must match at C-terminus (leader_core)
  or N-terminus (core_leader) depending on --core-orientation.
"""

from __future__ import annotations

import argparse
import hashlib
import logging
import re
from pathlib import Path
from statistics import median
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lasso_workbench.pipeline.orf_extraction import chunk_orfs
from lasso_workbench.utils.sequence_io import iter_genbank_records

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lab_csv", type=Path, required=True)
    parser.add_argument("--genomes_dir", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--window_nt", type=int, default=20000)
    parser.add_argument("--min_aa", type=int, default=20)
    parser.add_argument("--max_aa", type=int, default=120)
    parser.add_argument("--core_suffix_slack", type=int, default=0)
    parser.add_argument(
        "--core-orientation",
        choices=["leader_core", "core_leader"],
        default="leader_core",
        help="Match core at C-terminus (leader_core) or N-terminus (core_leader)",
    )
    parser.add_argument("--ambig_distance_nt", type=int, default=500)
    parser.add_argument("--limit_cases", type=int, default=0)
    parser.add_argument("--log_level", default="INFO")
    return parser.parse_args()


def normalize_accession(value: str) -> str:
    raw = str(value or "").strip()
    if not raw or raw.lower() == "unknown":
        return ""
    # Remove trailing comments and suffix letters (e.g. WP_...1c)
    raw = raw.split()[0].strip()
    if "." in raw:
        base, suffix = raw.rsplit(".", 1)
        if suffix.isdigit():
            raw = base
    return raw


def normalize_il(seq: str) -> str:
    return seq.upper().replace("I", "L")


def sanitize_core(seq: str) -> str:
    core = re.sub(r"[^A-Za-z]", "", str(seq or "")).upper()
    return core


def core_group_id(core: str) -> str:
    key = normalize_il(core)
    digest = hashlib.sha1(key.encode("utf-8")).hexdigest()[:10]
    return f"core_{digest}"


def distance_to_interval(start: int, end: int, anchor_start: int, anchor_end: int) -> int:
    if end <= anchor_start:
        return anchor_start - end
    if start >= anchor_end:
        return start - anchor_end
    return 0


def match_core_suffix(seq: str, core: str, tail_slack: int) -> bool:
    """Match core at the C-terminus (allowing a small trailing slack)."""
    if not core:
        return False
    seq_norm = normalize_il(seq)
    core_norm = normalize_il(core)
    n = len(core_norm)
    if n == 0 or len(seq_norm) < n:
        return False
    if tail_slack <= 0:
        return seq_norm.endswith(core_norm)
    tail = seq_norm[-(n + tail_slack):]
    return core_norm in tail


def match_core_prefix(seq: str, core: str, head_slack: int) -> bool:
    """Match core at the N-terminus (allowing a small leading slack)."""
    if not core:
        return False
    seq_norm = normalize_il(seq)
    core_norm = normalize_il(core)
    n = len(core_norm)
    if n == 0 or len(seq_norm) < n:
        return False
    if head_slack <= 0:
        return seq_norm.startswith(core_norm)
    head = seq_norm[: n + head_slack]
    return core_norm in head


def build_peptidase_index(
    genomes_dir: Path, accessions: Iterable[str]
) -> Dict[str, List[dict]]:
    targets = {acc for acc in accessions if acc}
    hits: Dict[str, List[PeptidaseLocus]] = {acc: [] for acc in targets}

    gbk_paths = list(genomes_dir.glob("**/*.gb*"))
    LOGGER.info("Scanning %d GBK files for peptidase accessions", len(gbk_paths))

    for gbk_path in gbk_paths:
        try:
            for record in iter_genbank_records(gbk_path):
                for feat in record.features:
                    if feat.type != "CDS":
                        continue
                    protein_id = feat.qualifiers.get("protein_id", [None])[0]
                    protein_id_norm = normalize_accession(protein_id)
                    if protein_id_norm not in targets:
                        continue
                    locus = {
                        "accession": protein_id_norm,
                        "gbk_path": gbk_path,
                        "record_id": str(record.id),
                        "start": int(feat.location.start),
                        "end": int(feat.location.end),
                        "strand": "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "?",
                        "product": feat.qualifiers.get("product", [None])[0],
                        "locus_tag": feat.qualifiers.get("locus_tag", [None])[0],
                    }
                    hits[protein_id_norm].append(locus)
        except Exception as exc:
            LOGGER.warning("Failed to parse %s: %s", gbk_path, exc)
            continue
    matched = sum(1 for v in hits.values() if v)
    LOGGER.info("Matched %d/%d accessions in genomes", matched, len(targets))
    return hits


def evaluate_locus(
    record,
    locus: dict,
    core_seq: str,
    window_nt: int,
    min_aa: int,
    max_aa: int,
    core_suffix_slack: int,
    core_orientation: str,
) -> dict:
    start = locus["start"]
    end = locus["end"]
    window_start = max(0, start - window_nt)
    window_end = min(len(record.seq), end + window_nt)
    window_seq = record.seq[window_start:window_end]

    orfs = list(chunk_orfs(window_seq, "+", window_start, window_end, min_aa, max_aa))
    rc_seq = window_seq.reverse_complement()
    orfs.extend(chunk_orfs(rc_seq, "-", window_start, window_end, min_aa, max_aa))

    match_a_orfs = []
    distances = []
    for orf in orfs:
        if core_orientation == "core_leader":
            match = match_core_prefix(orf.protein, core_seq, head_slack=core_suffix_slack)
        else:
            match = match_core_suffix(orf.protein, core_seq, tail_slack=core_suffix_slack)
        if match:
            match_a_orfs.append(orf)
            distances.append(distance_to_interval(orf.genomic_start, orf.genomic_end, start, end))

    match_count = len(match_a_orfs)
    med_dist = float(median(distances)) if distances else float("inf")
    locus_span = 0
    if match_a_orfs:
        starts = [orf.genomic_start for orf in match_a_orfs]
        ends = [orf.genomic_end for orf in match_a_orfs]
        locus_span = max(ends) - min(starts)

    return {
        "locus": locus,
        "orfs": orfs,
        "match_a_orfs": match_a_orfs,
        "match_count": match_count,
        "median_distance": med_dist,
        "locus_span": locus_span,
        "window_start": window_start,
        "window_end": window_end,
    }

def rank_loci(
    entries: List[dict],
    ambig_distance: int,
) -> Tuple[List[dict], List[int]]:
    ranked = sorted(
        entries,
        key=lambda item: (
            -item["match_count"],
            item["median_distance"],
            item["locus_span"],
            str(item["locus"]["gbk_path"]),
            item["locus"]["record_id"],
            item["locus"]["start"],
            item["locus"]["end"],
        ),
    )
    if not ranked:
        return ranked, []

    best_result = ranked[0]

    ambiguous_indices: List[int] = []
    for idx, result in enumerate(ranked[1:], start=1):
        if result["match_count"] != best_result["match_count"]:
            continue
        if abs(result["median_distance"] - best_result["median_distance"]) <= ambig_distance:
            ambiguous_indices.append(idx)

    if ambiguous_indices:
        ambiguous_indices.insert(0, 0)
    return ranked, ambiguous_indices


def write_fasta(records: Iterable[Tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    df = pd.read_csv(args.lab_csv, sep=";", encoding="utf-8-sig")
    if args.limit_cases > 0:
        df = df.head(args.limit_cases)

    df["peptidase_accession"] = df["Peptidase accessione"].map(normalize_accession)
    df["core_aa"] = df["Core peptide sequence"].map(sanitize_core)

    accessions = df["peptidase_accession"].dropna().unique().tolist()
    hits = build_peptidase_index(args.genomes_dir, accessions)

    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    case_rows: List[Dict[str, object]] = []
    locus_rows: List[Dict[str, object]] = []
    candidate_rows: List[Dict[str, object]] = []
    fasta_records: List[Tuple[str, str]] = []
    seen_candidates: set[Tuple[str, str]] = set()

    stats = {
        "total_cases": 0,
        "mapped_peptidase": 0,
        "cases_with_match_A": 0,
        "cases_unambiguous": 0,
        "cases_ambiguous": 0,
        "cases_no_core_match": 0,
        "loci_scanned": 0,
        "loci_included": 0,
        "ambiguous_loci_included": 0,
        "core_orientation": args.core_orientation,
    }

    record_cache: Dict[Tuple[str, str], Optional[object]] = {}

    for idx, row in df.iterrows():
        row_idx = idx + 1
        case_id = f"lab_{row_idx:04d}"
        name = str(row.get("Lasso peptide", f"case_{row_idx}"))
        core_seq = row.get("core_aa", "")
        accession = row.get("peptidase_accession", "")
        stats["total_cases"] += 1

        if not accession or not core_seq:
            continue

        loci = hits.get(accession, [])
        if not loci:
            continue

        stats["mapped_peptidase"] += 1

        locus_entries: List[dict] = []
        for locus in loci:
            record_key = (str(locus["gbk_path"]), locus["record_id"])
            if record_key not in record_cache:
                record = None
                for rec in iter_genbank_records(locus["gbk_path"]):
                    if str(rec.id) == locus["record_id"]:
                        record = rec
                        break
                record_cache[record_key] = record
            record = record_cache[record_key]
            if record is None:
                continue
            locus_result = evaluate_locus(
                record,
                locus,
                core_seq,
                args.window_nt,
                args.min_aa,
                args.max_aa,
                args.core_suffix_slack,
                args.core_orientation,
            )
            locus_entries.append(locus_result)

        if not locus_entries:
            continue

        stats["loci_scanned"] += len(locus_entries)

        case_has_match_a = any(entry["match_count"] > 0 for entry in locus_entries)

        if case_has_match_a:
            stats["cases_with_match_A"] += 1

        if not case_has_match_a:
            stats["cases_no_core_match"] += 1
            case_rows.append(
                {
                    "case_id": case_id,
                    "peptide_name": name,
                    "core_aa": core_seq,
                    "core_group": core_group_id(core_seq),
                    "peptidase_accession": accession,
                    "loci_total": len(locus_entries),
                    "loci_with_match_A": 0,
                    "case_ambiguous": False,
                }
            )
            continue

        ranked_entries, ambiguous_indices = rank_loci(locus_entries, args.ambig_distance_nt)
        case_ambiguous = len(ambiguous_indices) > 1
        if case_ambiguous:
            stats["cases_ambiguous"] += 1
        else:
            stats["cases_unambiguous"] += 1

        loci_with_match_a = sum(1 for entry in locus_entries if entry["match_count"] > 0)

        case_rows.append(
            {
                "case_id": case_id,
                "peptide_name": name,
                "core_aa": core_seq,
                "core_group": core_group_id(core_seq),
                "peptidase_accession": accession,
                "loci_total": len(locus_entries),
                "loci_with_match_A": loci_with_match_a,
                "case_ambiguous": case_ambiguous,
            }
        )

        core_group = core_group_id(core_seq)
        for rank_idx, locus_result in enumerate(ranked_entries, start=1):
            if case_ambiguous and (rank_idx - 1) not in ambiguous_indices:
                continue
            if not case_ambiguous and rank_idx != 1:
                continue

            locus = locus_result["locus"]
            window_start = locus_result["window_start"]
            window_end = locus_result["window_end"]

            locus_id = f"{case_id}_locus_{rank_idx:02d}"
            locus_ambiguous = case_ambiguous and (rank_idx - 1) in ambiguous_indices

            locus_rows.append(
                {
                    "case_id": case_id,
                    "locus_id": locus_id,
                    "locus_rank": rank_idx,
                    "case_ambiguous": case_ambiguous,
                    "locus_ambiguous": locus_ambiguous,
                    "peptide_name": name,
                    "core_aa": core_seq,
                    "core_group": core_group,
                    "core_orientation": args.core_orientation,
                    "peptidase_accession": accession,
                    "assembly_accession": locus["gbk_path"].parent.name,
                    "contig": locus["record_id"],
                    "peptidase_start": locus["start"],
                    "peptidase_end": locus["end"],
                    "peptidase_strand": locus["strand"],
                    "window_start": window_start,
                    "window_end": window_end,
                    "gbk_path": str(locus["gbk_path"]),
                    "matches_A": locus_result["match_count"],
                    "candidate_count": len(locus_result["orfs"]),
                    "median_distance_nt": locus_result["median_distance"],
                    "locus_span_nt": locus_result["locus_span"],
                }
            )

            stats["loci_included"] += 1
            if locus_ambiguous:
                stats["ambiguous_loci_included"] += 1

            for orf_idx, orf in enumerate(locus_result["match_a_orfs"], start=1):
                cand_id = f"{core_group}|case={case_id}|locus={locus_id}|pep={accession}|orf={orf_idx:04d}"
                seq = orf.protein
                key = (core_group, seq)
                if key not in seen_candidates:
                    fasta_records.append((cand_id, seq))
                    seen_candidates.add(key)

                candidate_rows.append(
                    {
                        "case_id": case_id,
                        "locus_id": locus_id,
                        "peptide_name": name,
                        "core_aa": core_seq,
                        "core_group": core_group,
                        "core_orientation": args.core_orientation,
                        "peptidase_accession": accession,
                        "assembly_accession": locus["gbk_path"].parent.name,
                        "contig": locus["record_id"],
                        "peptidase_start": locus["start"],
                        "peptidase_end": locus["end"],
                        "peptidase_strand": locus["strand"],
                        "candidate_id": cand_id,
                        "genomic_start": int(orf.genomic_start),
                        "genomic_end": int(orf.genomic_end),
                        "strand": orf.strand,
                        "aa_length": int(orf.aa_len),
                        "dist_peptidase_nt": distance_to_interval(
                            orf.genomic_start, orf.genomic_end, locus["start"], locus["end"]
                        ),
                        "protein_sequence": seq,
                    }
                )

    suffix = "" if args.core_orientation == "leader_core" else f"_{args.core_orientation}"
    write_fasta(fasta_records, out_dir / f"precursor_proteins_multi_strict{suffix}.faa")
    pd.DataFrame(locus_rows).to_csv(out_dir / f"lab_core_loci_multi_strict{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(case_rows).to_csv(out_dir / f"lab_core_cases_summary{suffix}.tsv", sep="\t", index=False)
    pd.DataFrame(candidate_rows).to_csv(out_dir / f"lab_core_candidates_multi_strict{suffix}.tsv", sep="\t", index=False)

    import json

    summary_path = out_dir / "lab_core_dataset_summary.json"
    summary_path.write_text(json.dumps(stats, indent=2) + "\n")

    LOGGER.info("Wrote %d strict precursors", len(fasta_records))
    LOGGER.info("Wrote %d locus rows", len(locus_rows))
    LOGGER.info("Wrote %d case summary rows", len(case_rows))
    LOGGER.info("Summary: %s", stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
