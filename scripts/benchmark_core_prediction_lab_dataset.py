#!/usr/bin/env python3
"""Benchmark core prediction on lab multi-candidate precursor dataset.

This uses the same CorePredictor + rules as the pipeline and evaluates
how often the predicted cores match the known core sequences.

Inputs:
  - precursor_proteins_multi_strict.faa (precursor sequences)
  - lab_core_candidates_multi_strict.tsv (candidate metadata + core_aa)

Outputs:
  - core_prediction_lab_summary.json
  - core_prediction_lab_cases.tsv
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Dict, List, Optional

import pandas as pd
from Bio import SeqIO

import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from lasso_workbench.core.prediction import CorePredictor, RuleEngine
from lasso_workbench.schemas.config import RuleSet
from lasso_workbench.services.rules_service import RulesService

LOGGER = logging.getLogger(__name__)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--precursor_faa", type=Path, default=None)
    parser.add_argument("--candidate_tsv", type=Path, required=True)
    parser.add_argument("--out_dir", type=Path, required=True)
    parser.add_argument("--top_n", type=int, default=5)
    parser.add_argument("--min_leader", type=int, default=10)
    parser.add_argument("--min_core", type=int, default=13)
    parser.add_argument("--limit_candidates", type=int, default=0)
    parser.add_argument(
        "--keep_duplicates",
        action="store_true",
        help="Keep duplicate protein_sequence/core_aa rows (default: dedupe).",
    )
    parser.add_argument("--rules_path", type=Path, default=None)
    parser.add_argument("--log_level", default="INFO")
    return parser.parse_args()


def normalize_il(seq: str) -> str:
    return str(seq).upper().replace("I", "L")


def matches_core_exact(seq: str, core: str) -> bool:
    if not core:
        return False
    seq_norm = normalize_il(seq)
    core_norm = normalize_il(core)
    if len(seq_norm) != len(core_norm):
        return False
    return seq_norm == core_norm


def compute_topk(ranks: List[Optional[int]], ks: List[int]) -> Dict[int, float]:
    if not ranks:
        return {k: 0.0 for k in ks}
    valid = [r for r in ranks if r is not None]
    return {k: sum(1 for r in valid if r <= k) / len(ranks) for k in ks}


def load_fasta_sequences(path: Path) -> Dict[str, str]:
    sequences: Dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        sequences[str(record.id)] = str(record.seq).upper()
    return sequences


def main() -> int:
    args = parse_args()
    logging.basicConfig(level=getattr(logging, args.log_level))

    if not args.candidate_tsv.exists():
        raise SystemExit(f"Missing candidate TSV: {args.candidate_tsv}")

    args.out_dir.mkdir(parents=True, exist_ok=True)

    candidates = pd.read_csv(args.candidate_tsv, sep="\t")
    if args.limit_candidates > 0:
        candidates = candidates.head(args.limit_candidates)

    if "protein_sequence" in candidates.columns:
        candidates["sequence"] = candidates["protein_sequence"].astype(str)
        if not args.keep_duplicates:
            candidates = candidates.drop_duplicates(subset=["protein_sequence", "core_aa"])
    elif args.precursor_faa is not None:
        if not args.precursor_faa.exists():
            raise SystemExit(f"Missing precursor FASTA: {args.precursor_faa}")
        seqs = load_fasta_sequences(args.precursor_faa)
        candidates["sequence"] = candidates["candidate_id"].map(seqs)
        missing_seq = candidates["sequence"].isna().sum()
        if missing_seq:
            LOGGER.warning("Missing sequences for %d candidates in FASTA", missing_seq)
            candidates = candidates.dropna(subset=["sequence"])
    else:
        raise SystemExit("No protein_sequence column and no --precursor_faa provided.")

    if candidates.empty:
        raise SystemExit("No candidates available for benchmarking.")

    if args.rules_path:
        rules = RuleSet.model_validate_json(args.rules_path.read_text())
        LOGGER.info("Loaded rules from %s", args.rules_path)
    else:
        rules_service = RulesService()
        rules = rules_service.get_rules()
        LOGGER.info("Loaded rules from RulesService (defaults + user overrides)")
    rule_engine = RuleEngine(rules)
    predictor = CorePredictor(rule_engine=rule_engine)

    variants = {
        "rules_only": None,
    }
    ks = [1, 3, 5]

    ranks_by_variant: Dict[str, List[Optional[int]]] = {v: [] for v in variants}
    best_rank_by_group: Dict[str, Dict[str, Optional[int]]] = {v: {} for v in variants}
    best_rank_by_case: Dict[str, Dict[str, Optional[int]]] = {v: {} for v in variants}

    summary = {
        "candidates_total": int(len(candidates)),
        "core_groups_total": int(candidates["core_group"].nunique()),
        "cases_total": int(candidates["case_id"].nunique()),
        "variants": {},
    }
    for variant in variants:
        summary["variants"][variant] = {
            "topk_per_candidate": {},
            "topk_best_per_core_group": {},
            "topk_best_per_case": {},
        }

    case_rows: List[Dict[str, object]] = []
    miss_rows: List[Dict[str, object]] = []

    for row in candidates.itertuples(index=False):
        candidate_id = str(getattr(row, "candidate_id"))
        core_aa = str(getattr(row, "core_aa", "")).strip().upper()
        if not core_aa:
            continue
        sequence = str(getattr(row, "sequence", "")).strip().upper()
        if not sequence:
            continue

        case_id = str(getattr(row, "case_id", ""))
        core_group = str(getattr(row, "core_group", ""))

        for variant in variants:
            prediction = predictor.predict(
                sequence,
                top_n=args.top_n,
                min_leader=args.min_leader,
                min_core=args.min_core,
            )

            predictions = prediction.predictions or []
            top1 = predictions[0] if predictions else None

            rank = None
            true_pred = None
            for idx, pred in enumerate(predictions, start=1):
                if matches_core_exact(pred.core, core_aa):
                    rank = idx
                    true_pred = pred
                    break
            ranks_by_variant[variant].append(rank)

            if core_group:
                existing = best_rank_by_group[variant].get(core_group)
                if existing is None or (rank is not None and rank < existing):
                    best_rank_by_group[variant][core_group] = rank

            if case_id:
                existing_case = best_rank_by_case[variant].get(case_id)
                if existing_case is None or (rank is not None and rank < existing_case):
                    best_rank_by_case[variant][case_id] = rank

            case_rows.append(
                {
                    "candidate_id": candidate_id,
                    "case_id": case_id,
                    "core_group": core_group,
                    "core_aa": core_aa,
                    "variant": variant,
                    "precursor_len": len(sequence),
                    "top1_core": top1.core if top1 else "",
                    "top1_score": top1.score if top1 else None,
                    "rank_exact": rank,
                }
            )

            if rank != 1:
                miss_rows.append(
                    {
                        "candidate_id": candidate_id,
                        "case_id": case_id,
                        "core_group": core_group,
                        "core_aa": core_aa,
                        "variant": variant,
                        "precursor_len": len(sequence),
                        "top1_core": top1.core if top1 else "",
                        "top1_score": top1.score if top1 else None,
                        "top1_reasons": ";".join(top1.reasons) if top1 and top1.reasons else "",
                        "true_rank": rank,
                        "true_score": true_pred.score if true_pred else None,
                        "true_reasons": ";".join(true_pred.reasons) if true_pred and true_pred.reasons else "",
                    }
                )

    for variant in variants:
        summary["variants"][variant]["topk_per_candidate"]["exact"] = compute_topk(
            ranks_by_variant[variant], ks
        )
        group_ranks = list(best_rank_by_group[variant].values())
        summary["variants"][variant]["topk_best_per_core_group"]["exact"] = compute_topk(
            group_ranks, ks
        )
        case_ranks = list(best_rank_by_case[variant].values())
        summary["variants"][variant]["topk_best_per_case"]["exact"] = compute_topk(
            case_ranks, ks
        )

    (args.out_dir / "core_prediction_lab_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    pd.DataFrame(case_rows).to_csv(
        args.out_dir / "core_prediction_lab_cases.tsv", sep="\t", index=False
    )
    if miss_rows:
        pd.DataFrame(miss_rows).to_csv(
            args.out_dir / "core_prediction_lab_misses.tsv", sep="\t", index=False
        )

    LOGGER.info("Wrote lab core prediction summary to %s", args.out_dir / "core_prediction_lab_summary.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
