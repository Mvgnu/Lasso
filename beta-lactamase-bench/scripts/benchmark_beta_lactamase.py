#!/usr/bin/env python3
"""Holdout benchmark for Ambler class A beta-lactamases (protein ID holdout)."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

import sys

ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lasso_workbench.pipeline.semantic_pipeline import (
    run_semantic_pipeline,
    results_to_dataframe,
    extract_orfs,
)
from lasso_workbench.schemas.pipeline import RankingConfig, BGCSegment
from lasso_workbench.utils.sequence_io import iter_genbank_records


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbk_dir", type=Path, required=True)
    ap.add_argument("--gbk_list", type=Path, default=None, help="Optional newline-separated GBK paths to use.")
    ap.add_argument("--validated_faa", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--model_name", type=str, required=True)
    ap.add_argument("--device", type=str, default="cpu")
    ap.add_argument("--limit_gbks", type=int, default=0)
    ap.add_argument("--limit_loci", type=int, default=0)
    ap.add_argument("--min_aa", type=int, default=200)
    ap.add_argument("--max_aa", type=int, default=400)
    ap.add_argument("--window_nt", type=int, default=20000)
    ap.add_argument(
        "--full-genome",
        action="store_true",
        help="Scan full record instead of a window around the protein_id locus.",
    )
    return ap.parse_args()


def read_fasta(path: Path) -> List[Tuple[str, str]]:
    records: List[Tuple[str, str]] = []
    name = None
    seq_parts: List[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq_parts)))
                name = line[1:].split()[0]
                seq_parts = []
            else:
                seq_parts.append(line)
    if name is not None:
        records.append((name, "".join(seq_parts)))
    return records


def write_fasta(records: Iterable[Tuple[str, str]], path: Path) -> None:
    with path.open("w") as handle:
        for name, seq in records:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def topk_success(rank: Optional[int], ks: List[int]) -> Dict[int, int]:
    out: Dict[int, int] = {}
    for k in ks:
        out[k] = int(rank is not None and rank <= k)
    return out


def normalize_accession(value: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    raw = raw.split()[0].strip()
    if "." in raw:
        base, suffix = raw.rsplit(".", 1)
        if suffix.isdigit():
            raw = base
    return raw


def parse_protein_accessions(records: Iterable[Tuple[str, str]]) -> Dict[str, List[str]]:
    mapping: Dict[str, List[str]] = {}
    for name, seq in records:
        acc = ""
        for part in name.split("|"):
            if part.startswith("prot="):
                acc = normalize_accession(part.split("=", 1)[1])
                break
        if not acc:
            acc = normalize_accession(name)
        if not acc:
            continue
        mapping.setdefault(acc, []).append(seq)
    return mapping


def extract_protein_accession(name: str) -> str:
    acc = ""
    for part in name.split("|"):
        if part.startswith("prot="):
            acc = normalize_accession(part.split("=", 1)[1])
            break
    if not acc:
        acc = normalize_accession(name)
    return acc


def collect_gbk_paths(root: Path) -> List[Path]:
    patterns = ("**/*.gb", "**/*.gbk", "**/*.gbff", "**/*.genbank")
    paths: List[Path] = []
    for pattern in patterns:
        paths.extend(root.glob(pattern))
    return sorted({p for p in paths if p.is_file()})


def build_protein_index(
    gbk_paths: List[Path], accessions: Iterable[str]
) -> Dict[str, List[dict]]:
    targets = {acc for acc in accessions if acc}
    hits: Dict[str, List[dict]] = {acc: [] for acc in targets}

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
                    hits[protein_id_norm].append(
                        {
                            "accession": protein_id_norm,
                            "gbk_path": gbk_path,
                            "record_id": str(record.id),
                            "start": int(feat.location.start),
                            "end": int(feat.location.end),
                            "strand": "+" if feat.location.strand == 1 else "-" if feat.location.strand == -1 else "?",
                        }
                    )
        except Exception:
            continue
    return hits


def build_segment(record, gbk_path: Path, index: int, window_start: int, window_end: int) -> BGCSegment:
    sequence = str(record.seq)
    return BGCSegment(
        record_id=str(record.id),
        index=index,
        sequence=sequence,
        annotated_cds=[],
        is_lasso=False,
        source_file=str(gbk_path),
        metadata={"region_ranges": [{"start": window_start, "end": window_end}]},
    )


def main() -> int:
    args = parse_args()

    if args.device not in {"auto", "cpu", "cuda", "mps"}:
        raise ValueError(f"Unsupported device: {args.device}")

    validated_records = read_fasta(args.validated_faa)
    accession_to_seqs = parse_protein_accessions(validated_records)
    if not accession_to_seqs:
        raise SystemExit("No protein accessions found in validated FASTA ids (expected prot=...).")

    if args.gbk_list:
        gbk_paths = [Path(line.strip()) for line in args.gbk_list.read_text().splitlines() if line.strip()]
    else:
        gbk_paths = collect_gbk_paths(args.gbk_dir)
    if args.limit_gbks > 0:
        gbk_paths = gbk_paths[: args.limit_gbks]
    if not gbk_paths:
        raise SystemExit(f"No GBK files found under {args.gbk_dir}")

    protein_hits = build_protein_index(gbk_paths, accession_to_seqs.keys())
    loci: List[dict] = []
    for acc, loci_list in protein_hits.items():
        loci.extend(loci_list)
    if args.limit_loci > 0:
        loci = loci[: args.limit_loci]

    args.out_dir.mkdir(parents=True, exist_ok=True)
    topk = [1, 5, 10, 50]
    per_bgc_rows: List[Dict[str, object]] = []
    candidate_counts: List[int] = []
    all_df_rows: List[pd.DataFrame] = []
    loci_with_matches = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir_path = Path(tmpdir)
        record_cache: Dict[Tuple[str, str], object] = {}
        for idx, locus in enumerate(loci, start=1):
            gbk = Path(locus["gbk_path"])
            record_id = locus["record_id"]
            cache_key = (str(gbk), record_id)
            if cache_key not in record_cache:
                rec = None
                for record in iter_genbank_records(gbk):
                    if str(record.id) == record_id:
                        rec = record
                        break
                record_cache[cache_key] = rec
            record = record_cache[cache_key]
            if record is None:
                continue

            if args.full_genome:
                window_start = 0
                window_end = len(record.seq)
            else:
                window_start = max(0, int(locus["start"]) - args.window_nt)
                window_end = min(len(record.seq), int(locus["end"]) + args.window_nt)
            segment = build_segment(record, gbk, idx, window_start, window_end)

            orfs, _ = extract_orfs([segment], min_aa=args.min_aa, max_aa=args.max_aa)
            holdout_seqs = {
                seq.upper() for seq in accession_to_seqs.get(locus["accession"], [])
            }
            if not holdout_seqs:
                continue
            matches: List[str] = []
            for orf in orfs:
                seq = orf.protein_sequence.upper()
                if seq in holdout_seqs:
                    matches.append(seq)

            if not matches:
                continue

            loci_with_matches += 1
            match_set = {s.upper() for s in matches}
            filtered = [
                (name, seq)
                for name, seq in validated_records
                if extract_protein_accession(name) != locus["accession"]
            ]
            if not filtered:
                continue

            holdout_faa = tmpdir_path / f"validated_holdout_{idx:05d}.faa"
            write_fasta(filtered, holdout_faa)

            results, df = run_semantic_pipeline(
                gbk_files=None,
                segments=[segment],
                validated_faa=holdout_faa,
                output_dir=args.out_dir,
                ranking_config=RankingConfig(),
                min_aa=args.min_aa,
                max_aa=args.max_aa,
                device=None if args.device == "auto" else args.device,
                model_name=args.model_name,
                progress=None,
            )
            df = results_to_dataframe(results)
            if df.empty:
                continue
            sorted_group = df.sort_values("best_similarity", ascending=False).reset_index(drop=True)
            sorted_group["exact_rank"] = sorted_group.index + 1
            holdout_mask = sorted_group["protein_sequence"].str.upper().isin(match_set)
            match_count = int(holdout_mask.sum())
            if match_count:
                rank = int(sorted_group.loc[holdout_mask, "exact_rank"].min())
            else:
                rank = None

            holdout_rows = sorted_group[holdout_mask].copy()
            if not holdout_rows.empty:
                holdout_rows["is_holdout"] = True
                holdout_rows["gbk_path"] = str(gbk)
                holdout_rows["window_start"] = window_start
                holdout_rows["window_end"] = window_end
                holdout_rows["protein_accession"] = locus["accession"]
                holdout_rows["rank"] = rank
                if "exact_rank" in holdout_rows.columns:
                    holdout_rows = holdout_rows.drop(columns=["exact_rank"])
                all_df_rows.append(holdout_rows)

            success = topk_success(rank, topk)
            gap = None
            if len(sorted_group) >= 2:
                gap = float(sorted_group.iloc[0].best_similarity - sorted_group.iloc[1].best_similarity)

            candidate_count = len(sorted_group)
            candidate_counts.append(candidate_count)
            per_bgc_rows.append(
                {
                    "record_id": record_id,
                    "rank": rank,
                    "top1": success[1],
                    "top5": success[5],
                    "top10": success[10],
                    "top50": success[50],
                    "top1_top2_gap": gap,
                    "candidate_count": candidate_count,
                    "holdout_match_count": match_count,
                }
            )

    per_bgc_df = pd.DataFrame(per_bgc_rows)
    per_bgc_df.to_csv(args.out_dir / "holdout_topk.tsv", sep="\t", index=False)

    holdout_df = pd.concat(all_df_rows, ignore_index=True) if all_df_rows else pd.DataFrame()
    if not holdout_df.empty:
        holdout_df.to_csv(args.out_dir / "holdout_candidates.tsv", sep="\t", index=False)

    summary = {
        "gbks_total": len(gbk_paths),
        "loci_total": len(loci),
        "loci_with_validated_protein": int(loci_with_matches),
        "loci_evaluated": int(per_bgc_df.shape[0]),
        "topk": {
            "top1": float(per_bgc_df["top1"].mean()) if not per_bgc_df.empty else 0.0,
            "top5": float(per_bgc_df["top5"].mean()) if not per_bgc_df.empty else 0.0,
            "top10": float(per_bgc_df["top10"].mean()) if not per_bgc_df.empty else 0.0,
            "top50": float(per_bgc_df["top50"].mean()) if not per_bgc_df.empty else 0.0,
        },
        "candidate_count": {
            "mean": float(sum(candidate_counts) / len(candidate_counts)) if candidate_counts else 0.0,
            "median": float(per_bgc_df["candidate_count"].median()) if not per_bgc_df.empty else 0.0,
            "min": int(per_bgc_df["candidate_count"].min()) if not per_bgc_df.empty else 0,
            "max": int(per_bgc_df["candidate_count"].max()) if not per_bgc_df.empty else 0,
        },
        "model_name": args.model_name,
        "device": args.device,
    }

    (args.out_dir / "holdout_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
