#!/usr/bin/env python3
"""
Fetch Ambler class A beta-lactamase protein FASTAs from NCBI using GenPept IDs.
License: MIT
Author: Magnus Ohle
"""
from __future__ import annotations

import argparse
import csv
import time
from pathlib import Path
from typing import Dict, List, Optional
from urllib.parse import urlencode
from urllib.request import urlopen


def normalize_accession(value: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    raw = raw.split()[0].strip()
    return raw.replace("\u00a0", "").strip()


def extract_accession(row: Dict[str, str]) -> str:
    acc = normalize_accession(row.get("genpept_id", ""))
    if acc:
        return acc
    link = row.get("genpept_link", "")
    if link:
        acc = normalize_accession(link.rstrip("/").split("/")[-1])
    return acc


def fetch_fasta_batch(accessions: List[str], email: str, tool: str, api_key: str) -> str:
    params = {
        "db": "protein",
        "id": ",".join(accessions),
        "rettype": "fasta",
        "retmode": "text",
    }
    if tool:
        params["tool"] = tool
    if email:
        params["email"] = email
    if api_key:
        params["api_key"] = api_key
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?" + urlencode(params)
    with urlopen(url) as response:  # nosec - NCBI fetch
        return response.read().decode("utf-8")


def parse_fasta_text(text: str) -> List[tuple[str, str]]:
    records: List[tuple[str, str]] = []
    header = ""
    seq_parts: List[str] = []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line)
    if header:
        records.append((header, "".join(seq_parts)))
    return records


def extract_header_accession(header: str) -> str:
    token = header.split()[0]
    acc = normalize_accession(token)
    if "." in acc:
        base, suffix = acc.rsplit(".", 1)
        if suffix.isdigit():
            acc = base
    return acc


def write_fasta(records: List[tuple[str, str]], output_path: Path) -> None:
    with output_path.open("w") as handle:
        for header, seq in records:
            handle.write(f">{header}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--table-tsv",
        type=Path,
        default=Path("beta-lactamase-bench/data/ambler_class_a_table.tsv"),
    )
    ap.add_argument(
        "--output-faa",
        type=Path,
        default=Path("beta-lactamase-bench/data/ambler_class_a_beta_lactamases.faa"),
    )
    ap.add_argument(
        "--output-tsv",
        type=Path,
        default=Path("beta-lactamase-bench/data/ambler_class_a_beta_lactamases.tsv"),
    )
    ap.add_argument("--email", type=str, default="magnus.ohle@student.uni-tuebingen.de")
    ap.add_argument("--tool", type=str, default="lasso_workbench")
    ap.add_argument("--sleep", type=float, default=0.5)
    ap.add_argument("--batch-size", type=int, default=200)
    ap.add_argument("--api-key", type=str, default="")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    if not args.table_tsv.exists():
        raise SystemExit(f"Missing table TSV: {args.table_tsv}")

    rows = []
    with args.table_tsv.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append(row)

    records: List[tuple[str, str]] = []
    output_rows: List[Dict[str, str]] = []
    seen = set()
    acc_to_row: Dict[str, Dict[str, str]] = {}
    accessions: List[str] = []

    for row in rows:
        accession = extract_accession(row)
        if not accession:
            continue
        if accession not in acc_to_row:
            acc_to_row[accession] = row
        if accession in seen:
            continue
        accessions.append(accession)
        seen.add(accession)

    if args.limit:
        accessions = accessions[: args.limit]

    batch_size = max(1, int(args.batch_size))
    for start in range(0, len(accessions), batch_size):
        batch = accessions[start : start + batch_size]
        print(f"Fetching {len(batch)} accessions ({start + 1}-{start + len(batch)} of {len(accessions)})")
        fasta_text = fetch_fasta_batch(batch, args.email, args.tool, args.api_key)
        if not fasta_text:
            continue
        for header_line, seq in parse_fasta_text(fasta_text):
            if not seq:
                continue
            accession = extract_header_accession(header_line)
            row = acc_to_row.get(accession, {})
            name = row.get("protein_name", "")
            header = (
                f"prot={accession}|locus={accession}|name={name}"
                f"|class={row.get('ambler_class','A')}|source=ambler_table"
            )
            records.append((header, seq))
            output_rows.append(
                {
                    "accession": accession,
                    "protein_name": name,
                    "sequence_length": str(len(seq)),
                    "header_line": header_line,
                    "genbank_id": row.get("genbank_id", ""),
                    "genbank_link": row.get("genbank_link", ""),
                    "subfamily": row.get("subfamily", ""),
                    "phenotype": row.get("phenotype", ""),
                    "natural_or_acquired": row.get("natural_or_acquired", ""),
                }
            )
        if args.sleep:
            time.sleep(args.sleep)

    if not records:
        raise SystemExit("No FASTA records fetched.")

    args.output_faa.parent.mkdir(parents=True, exist_ok=True)
    write_fasta(records, args.output_faa)
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(output_rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"Wrote {len(records)} sequences to {args.output_faa}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
