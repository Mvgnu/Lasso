#!/usr/bin/env python3
"""Per-GBK holdout benchmark for pipeline ranking (leave-one-out)."""

from __future__ import annotations

import argparse
import json
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd

import sys

_ROOT = Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from lasso_workbench.pipeline.semantic_pipeline import (
    run_semantic_pipeline,
    parse_gbk_file,
    results_to_dataframe,
    extract_orfs,
)
from lasso_workbench.schemas.pipeline import RankingConfig


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--gbk_dir", type=Path, required=True)
    ap.add_argument("--gbk_list", type=Path, default=None, help="Optional newline-separated GBK paths to use.")
    ap.add_argument("--validated_faa", type=Path, required=True)
    ap.add_argument("--out_dir", type=Path, required=True)
    ap.add_argument("--model_name", type=str, required=True)
    ap.add_argument("--device", type=str, default="cpu")
    ap.add_argument("--limit_gbks", type=int, default=0)
    ap.add_argument("--min_aa", type=int, default=20)
    ap.add_argument("--max_aa", type=int, default=120)
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


def main() -> int:
    args = parse_args()

    if args.device not in {"auto", "cpu", "cuda", "mps"}:
        raise ValueError(f"Unsupported device: {args.device}")

    if args.gbk_list:
        gbk_paths = [Path(line.strip()) for line in args.gbk_list.read_text().splitlines() if line.strip()]
    else:
        gbk_paths = sorted(args.gbk_dir.glob("**/*.gb*"))
    if args.limit_gbks > 0:
        gbk_paths = gbk_paths[: args.limit_gbks]

    validated_records = read_fasta(args.validated_faa)
    validated_seqs = {seq.upper() for _, seq in validated_records}

    args.out_dir.mkdir(parents=True, exist_ok=True)
    topk = [1, 5, 10]
    per_bgc_rows: List[Dict[str, object]] = []
    all_df_rows: List[pd.DataFrame] = []
    gbks_with_matches = 0

    with tempfile.TemporaryDirectory() as tmpdir:
        holdout_faa = Path(tmpdir) / "validated_holdout.faa"

        for idx, gbk in enumerate(gbk_paths, start=1):
            segment = parse_gbk_file(gbk, idx)
            if not segment:
                continue
            orfs, _ = extract_orfs([segment], min_aa=args.min_aa, max_aa=args.max_aa)
            matches: List[str] = []
            for orf in orfs:
                seq = orf.protein_sequence.upper()
                if seq in validated_seqs:
                    matches.append(seq)

            if not matches:
                continue

            gbks_with_matches += 1
            match_set = {s.upper() for s in matches}
            filtered = [(name, seq) for name, seq in validated_records if seq.upper() not in match_set]
            if not filtered:
                continue

            write_fasta(filtered, holdout_faa)

            results, df = run_semantic_pipeline(
                gbk_files=[gbk],
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
            all_df_rows.append(df)

            record_id = segment.record_id
            sorted_group = df.sort_values("best_similarity", ascending=False)
            rank = None
            for idx_row, row in enumerate(sorted_group.itertuples(index=False), start=1):
                if str(row.protein_sequence).upper() in match_set:
                    rank = idx_row
                    break

            success = topk_success(rank, topk)
            gap = None
            if len(sorted_group) >= 2:
                gap = float(sorted_group.iloc[0].best_similarity - sorted_group.iloc[1].best_similarity)

            per_bgc_rows.append(
                {
                    "record_id": record_id,
                    "rank": rank,
                    "top1": success[1],
                    "top5": success[5],
                    "top10": success[10],
                    "top1_top2_gap": gap,
                    "candidate_count": len(sorted_group),
                }
            )

    per_bgc_df = pd.DataFrame(per_bgc_rows)
    per_bgc_df.to_csv(args.out_dir / "holdout_topk.tsv", sep="\t", index=False)

    if all_df_rows:
        pd.concat(all_df_rows, ignore_index=True).to_csv(
            args.out_dir / "holdout_candidates.tsv", sep="\t", index=False
        )

    summary = {
        "gbks_total": len(gbk_paths),
        "gbks_with_validated_precursor": int(gbks_with_matches),
        "gbks_evaluated": int(per_bgc_df.shape[0]),
        "topk": {
            "top1": float(per_bgc_df["top1"].mean()) if not per_bgc_df.empty else 0.0,
            "top5": float(per_bgc_df["top5"].mean()) if not per_bgc_df.empty else 0.0,
            "top10": float(per_bgc_df["top10"].mean()) if not per_bgc_df.empty else 0.0,
        },
        "model_name": args.model_name,
        "device": args.device,
    }

    (args.out_dir / "holdout_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
