#!/usr/bin/env python3
"""Score provided precursor/core pairs with the rule engine and report core rank."""
from __future__ import annotations

import argparse
import csv
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Tuple

from lasso_workbench.core.prediction import CorePredictor


@dataclass
class Record:
    name: str
    precursor: str
    core: str


def _clean_seq(text: str) -> str:
    return "".join(re.findall(r"[A-Za-z]", text)).upper()


def _parse_file(path: Path) -> List[Record]:
    records: List[Record] = []
    name = ""
    precursor = ""
    core = ""
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and precursor and core:
                records.append(Record(name=name, precursor=precursor, core=core))
            name = line[1:].strip()
            precursor = ""
            core = ""
            continue
        if line.lower().startswith("precursor"):
            line = line.split(":", 1)[-1].strip() if ":" in line else line[len("precursor"):].strip()
            precursor = _clean_seq(line)
            continue
        if line.lower().startswith("core"):
            line = line.split(":", 1)[-1].strip() if ":" in line else line[len("core"):].strip()
            core = _clean_seq(line)
            continue
    if name and precursor and core:
        records.append(Record(name=name, precursor=precursor, core=core))
    return records


def _allowed_positions(seq_len: int, min_leader: int, min_core: int, orientation: str) -> List[int]:
    if orientation == "core_leader":
        start = min_core
        end = seq_len - min_leader
    else:
        start = min_leader
        end = seq_len - min_core
    if end < start:
        return []
    return list(range(start, end + 1))


def _score_positions(
    predictor: CorePredictor,
    sequence: str,
    min_leader: int,
    min_core: int,
    orientation: str,
) -> List[Tuple[int, float]]:
    positions = _allowed_positions(len(sequence), min_leader, min_core, orientation)
    scores: List[Tuple[int, float]] = []
    for pos in positions:
        score, _ = predictor.rules.score_cleavage_site(sequence, pos, orientation=orientation)
        scores.append((pos, float(score)))
    scores.sort(key=lambda x: (-x[1], x[0]))
    return scores


def _find_core_positions(sequence: str, core: str, orientation: str) -> List[int]:
    if not core:
        return []
    if orientation == "core_leader":
        if sequence.startswith(core):
            return [len(core)]
        return []
    positions: List[int] = []
    start = 0
    while True:
        idx = sequence.find(core, start)
        if idx == -1:
            break
        positions.append(idx)
        start = idx + 1
    return positions


def main() -> int:
    parser = argparse.ArgumentParser(description="Score newlassos precursors with rule engine.")
    parser.add_argument("--input", type=Path, default=Path("data/lab_dataset/newlassos.txt"))
    parser.add_argument("--output", type=Path, default=Path("results/newlassos_rule_scores.tsv"))
    parser.add_argument("--min-leader", type=int, default=10)
    parser.add_argument("--min-core", type=int, default=None)
    parser.add_argument("--top-n", type=int, default=3)
    parser.add_argument("--invert-names", type=str, default="", help="Comma-separated substrings to score in core_leader orientation")
    args = parser.parse_args()

    records = _parse_file(args.input)
    if not records:
        raise SystemExit("No records found in input file.")

    predictor = CorePredictor()
    if args.min_core is None:
        rule = predictor.rules._get_rule("core_length")
        if rule and rule.enabled:
            args.min_core = rule.parameters.min
        else:
            args.min_core = 10

    rows = []
    invert_tokens = [s.strip().lower() for s in args.invert_names.split(",") if s.strip()]

    for rec in records:
        orientation = "leader_core"
        if invert_tokens and any(tok in rec.name.lower() for tok in invert_tokens):
            orientation = "core_leader"

        preds = predictor.predict(rec.precursor, top_n=args.top_n, allow_inverted=orientation == "core_leader")
        if orientation == "core_leader":
            best = next((p for p in preds.predictions if p.orientation == "core_leader"), None)
        else:
            best = preds.best_prediction

        scores = _score_positions(predictor, rec.precursor, args.min_leader, args.min_core, orientation)
        pos_to_rank = {pos: idx + 1 for idx, (pos, _) in enumerate(scores)}
        pos_to_score = {pos: score for pos, score in scores}

        core_positions = _find_core_positions(rec.precursor, rec.core, orientation)
        chosen_pos: Optional[int] = None
        chosen_score: Optional[float] = None
        chosen_rank: Optional[int] = None
        if core_positions:
            # pick the core start with the best rule score
            best_pos = max(core_positions, key=lambda p: pos_to_score.get(p, float("-inf")))
            chosen_pos = best_pos
            chosen_score = pos_to_score.get(best_pos)
            chosen_rank = pos_to_rank.get(best_pos)

        rows.append({
            "name": rec.name,
            "precursor": rec.precursor,
            "core_given": rec.core,
            "orientation": orientation,
            "core_positions": ",".join(str(p) for p in core_positions) if core_positions else "",
            "core_score": f"{chosen_score:.2f}" if chosen_score is not None else "",
            "core_rank": chosen_rank or "",
            "core_rank_total": len(scores),
            "best_pred_core": best.core if best else "",
            "best_pred_pos": best.cleavage_site if best else "",
            "best_pred_score": f"{best.score:.2f}" if best else "",
            "best_pred_confidence": best.confidence if best else "",
        })

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
