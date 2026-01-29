#!/usr/bin/env python3
"""Build a SQLite index of 20-120 aa ORFs from antiSMASH lasso regions."""
from __future__ import annotations

import argparse
import logging
from pathlib import Path

from lasso_workbench.altframe import build_lasso_orf_index


def main() -> int:
    parser = argparse.ArgumentParser(description="Build SQLite ORF index for lasso-annotated regions.")
    parser.add_argument(
        "--gbk-dir",
        type=Path,
        default=Path("data/antismash_lasso/gbk"),
        help="Directory with antiSMASH GBK files.",
    )
    parser.add_argument(
        "--db-path",
        type=Path,
        default=Path("results/lasso_orf_index.sqlite"),
        help="SQLite database path to create/update.",
    )
    parser.add_argument("--min-aa", type=int, default=20, help="Minimum amino acid length.")
    parser.add_argument("--max-aa", type=int, default=120, help="Maximum amino acid length.")
    parser.add_argument(
        "--reset",
        action="store_true",
        help="Drop existing tables before indexing.",
    )
    parser.add_argument(
        "--batch-size",
        type=int,
        default=5000,
        help="Batch size for DB inserts.",
    )
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    stats = build_lasso_orf_index(
        gbk_dir=args.gbk_dir,
        db_path=args.db_path,
        min_aa=args.min_aa,
        max_aa=args.max_aa,
        reset=args.reset,
        batch_size=args.batch_size,
    )

    logging.info("Index complete: %s", stats)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
