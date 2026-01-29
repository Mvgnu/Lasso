#!/usr/bin/env python3
"""
Inject locus= tokens into FASTA headers using prot= (or first token).
License: MIT
Author: Magnus Ohle
"""
from __future__ import annotations

import argparse
from pathlib import Path


def derive_locus(header: str) -> str:
    for part in header.split("|"):
        if part.startswith("prot="):
            return part.split("=", 1)[1]
    return header.split()[0]


def patch_headers(input_path: Path, output_path: Path) -> None:
    with input_path.open() as handle:
        lines = handle.read().splitlines()

    out_lines: list[str] = []
    for line in lines:
        if not line.startswith(">"):
            out_lines.append(line)
            continue
        header = line[1:].strip()
        if "locus=" in header:
            out_lines.append(">" + header)
            continue
        locus = derive_locus(header)
        if "|" in header:
            header = header.replace("|", f"|locus={locus}|", 1)
        else:
            header = f"{header}|locus={locus}"
        out_lines.append(">" + header)

    output_path.write_text("\n".join(out_lines) + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--input-faa",
        type=Path,
        default=Path("beta-lactamase-bench/data/ambler_class_a_beta_lactamases.faa"),
    )
    ap.add_argument("--output-faa", type=Path, default=None)
    args = ap.parse_args()

    output_path = args.output_faa or args.input_faa
    if not args.input_faa.exists():
        raise SystemExit(f"Missing input FASTA: {args.input_faa}")
    patch_headers(args.input_faa, output_path)
    print(f"Wrote {output_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
