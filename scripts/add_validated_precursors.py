#!/usr/bin/env python3
"""
Add curated precursor sequences to a validated dataset with consistent IDs.
License: MIT
Author: Magnus Ohle

Creates/updates:
- data/precursors/precursor_proteins_curated.faa
- data/precursors/precursor_proteins_curated.tsv
Optionally writes a merged validated FASTA combining the generated multi_strict set
with curated additions.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO


@dataclass
class Entry:
    name: str
    sequence: str
    core: str = ""
    orientation: str = "leader_core"


def _clean_seq(text: str) -> str:
    return "".join(re.findall(r"[A-Za-z]", text)).upper()


def _slug(text: str) -> str:
    text = text.strip().lower()
    text = re.sub(r"\s+", "_", text)
    text = re.sub(r"[^a-z0-9_.-]", "", text)
    return text[:60] if text else "unnamed"


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:10]


def _parse_newlassos(path: Path) -> List[Entry]:
    entries: List[Entry] = []
    name = ""
    precursor = ""
    core = ""
    for raw in path.read_text().splitlines():
        line = raw.strip()
        if not line:
            continue
        if line.startswith(">"):
            if name and precursor:
                entries.append(Entry(name=name, sequence=precursor, core=core))
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
    if name and precursor:
        entries.append(Entry(name=name, sequence=precursor, core=core))
    return entries


def _parse_fasta(path: Path, orientation: str) -> List[Entry]:
    entries: List[Entry] = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = _clean_seq(str(record.seq))
        if not seq:
            continue
        desc = record.description or record.id
        name = desc
        core = ""
        match = re.search(r"core[:=]([A-Za-z]+)", desc)
        if match:
            core = _clean_seq(match.group(1))
        entries.append(Entry(name=name, sequence=seq, core=core, orientation=orientation))
    return entries


def _first_nonempty(row: dict, keys: List[str]) -> str:
    for key in keys:
        value = row.get(key)
        if value is None:
            continue
        text = str(value).strip()
        if text and text.lower() != "nan":
            return text
    return ""


def _parse_table(path: Path, orientation: str) -> List[Entry]:
    df = pd.read_csv(path, sep="\t" if path.suffix.lower() == ".tsv" else ",")
    if df.empty:
        return []
    lower_map = {col: col.strip().lower() for col in df.columns}
    df = df.rename(columns=lower_map)

    entries: List[Entry] = []
    for _, row in df.iterrows():
        row_dict = row.to_dict()
        name = _first_nonempty(row_dict, ["name", "id", "label", "peptide_name"])
        seq = _first_nonempty(row_dict, ["sequence", "precursor", "full_precursor", "protein_sequence"])
        core = _first_nonempty(row_dict, ["core", "core_aa", "core_sequence"])
        orient = _first_nonempty(row_dict, ["orientation", "core_orientation"]) or orientation
        seq = _clean_seq(seq)
        core = _clean_seq(core)
        if not seq:
            continue
        entries.append(Entry(name=name or f"seq_{_seq_hash(seq)}", sequence=seq, core=core, orientation=orient))
    return entries


def _read_existing_faa(path: Path) -> Dict[str, str]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    seqs: Dict[str, str] = {}
    for record in SeqIO.parse(str(path), "fasta"):
        seq = _clean_seq(str(record.seq))
        if not seq:
            continue
        seqs[record.id] = seq
    return seqs


def _read_existing_tsv(path: Path) -> List[dict]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _write_faa(path: Path, entries: List[Tuple[str, str]]) -> None:
    with path.open("w") as handle:
        for name, seq in entries:
            handle.write(f">{name}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")

def _write_verified_tsv(
    path: Path,
    multi_rows: List[dict],
    curated_rows: List[dict],
    multi_strict_source: Path,
) -> None:
    fieldnames = [
        "id",
        "name",
        "sequence",
        "length",
        "core",
        "orientation",
        "locus",
        "source",
        "added_at",
        "input",
    ]
    seen: set[str] = set()
    rows: List[dict] = []

    for row in multi_rows:
        seq = _clean_seq(row.get("protein_sequence", ""))
        if not seq or seq in seen:
            continue
        rows.append(
            {
                "id": row.get("candidate_id") or f"multi_{_seq_hash(seq)}",
                "name": row.get("peptide_name") or row.get("candidate_id") or "",
                "sequence": seq,
                "length": row.get("aa_length") or len(seq),
                "core": _clean_seq(row.get("core_aa", "")),
                "orientation": "leader_core",
                "locus": row.get("locus_id") or "",
                "source": "multi_strict",
                "added_at": "",
                "input": str(multi_strict_source),
            }
        )
        seen.add(seq)

    for row in curated_rows:
        seq = _clean_seq(row.get("sequence", ""))
        if not seq or seq in seen:
            continue
        rows.append(
            {
                "id": row.get("id") or f"curated_{_seq_hash(seq)}",
                "name": row.get("name") or "",
                "sequence": seq,
                "length": row.get("length") or len(seq),
                "core": _clean_seq(row.get("core", "")),
                "orientation": row.get("orientation") or "leader_core",
                "locus": row.get("locus") or "",
                "source": row.get("source") or "curated",
                "added_at": row.get("added_at") or "",
                "input": row.get("input") or "",
            }
        )
        seen.add(seq)

    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--format", choices=["auto", "fasta", "newlassos", "tsv", "csv"], default="auto")
    parser.add_argument("--output-dir", type=Path, default=Path("data/precursors"))
    parser.add_argument("--curated-faa", type=str, default="precursor_proteins_curated.faa")
    parser.add_argument("--curated-tsv", type=str, default="precursor_proteins_curated.tsv")
    parser.add_argument("--merge", action="store_true", default=True)
    parser.add_argument("--no-merge", action="store_false", dest="merge")
    parser.add_argument("--merge-output", type=str, default="precursor_proteins_verified.faa")
    parser.add_argument("--merge-tsv-output", type=str, default="precursor_proteins_verified.tsv")
    parser.add_argument("--source", type=str, default="curated")
    parser.add_argument("--orientation", choices=["leader_core", "core_leader"], default="leader_core")
    args = parser.parse_args()

    fmt = args.format
    if fmt == "auto":
        suffix = args.input.suffix.lower()
        if suffix in {".faa", ".fasta", ".fa"}:
            fmt = "fasta"
        elif suffix in {".tsv", ".csv"}:
            fmt = suffix.lstrip(".")
        else:
            fmt = "fasta"
    if fmt == "newlassos":
        new_entries = _parse_newlassos(args.input)
    elif fmt in {"tsv", "csv"}:
        new_entries = _parse_table(args.input, args.orientation)
    else:
        new_entries = _parse_fasta(args.input, args.orientation)
    if not new_entries:
        raise SystemExit("No sequences parsed from input.")

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    curated_faa = out_dir / args.curated_faa
    curated_tsv = out_dir / args.curated_tsv

    existing_faa = _read_existing_faa(curated_faa)
    existing_tsv = _read_existing_tsv(curated_tsv)
    existing_ids = {row.get("id", "") for row in existing_tsv}
    existing_seqs = {seq for seq in existing_faa.values()}

    added = 0
    now = datetime.now(timezone.utc).isoformat()
    updated_rows = existing_tsv[:]
    updated_entries = list(existing_faa.items())

    for entry in new_entries:
        seq = entry.sequence
        if not seq:
            continue
        if seq in existing_seqs:
            continue
        seq_hash = _seq_hash(seq)
        locus = f"curated_{seq_hash}"
        name_slug = _slug(entry.name)
        rec_id = f"curated_{seq_hash}|locus={locus}|name={name_slug}|source={args.source}"
        if rec_id in existing_ids:
            continue
        updated_entries.append((rec_id, seq))
        updated_rows.append({
            "id": rec_id,
            "name": entry.name,
            "sequence": seq,
            "length": len(seq),
            "core": entry.core or "",
            "orientation": entry.orientation,
            "locus": locus,
            "source": args.source,
            "added_at": now,
            "input": str(args.input),
        })
        existing_seqs.add(seq)
        existing_ids.add(rec_id)
        added += 1

    _write_faa(curated_faa, updated_entries)
    with curated_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["id", "name", "sequence", "length", "core", "orientation", "locus", "source", "added_at", "input"],
            delimiter="\t",
        )
        writer.writeheader()
        writer.writerows(updated_rows)

    if args.merge:
        generated = out_dir / "precursor_proteins_multi_strict.faa"
        multi_strict_tsv = out_dir / "lab_core_candidates_multi_strict.tsv"
        merged_entries: List[Tuple[str, str]] = []
        seen = set()
        for path in (generated, curated_faa):
            if not path.exists():
                continue
            for record in SeqIO.parse(str(path), "fasta"):
                seq = _clean_seq(str(record.seq))
                if not seq or seq in seen:
                    continue
                merged_entries.append((record.id, seq))
                seen.add(seq)
        merged_path = out_dir / args.merge_output
        _write_faa(merged_path, merged_entries)
        if multi_strict_tsv.exists():
            multi_rows = _read_existing_tsv(multi_strict_tsv)
        else:
            multi_rows = []
        verified_tsv_path = out_dir / args.merge_tsv_output
        _write_verified_tsv(verified_tsv_path, multi_rows, updated_rows, multi_strict_tsv)

    print(f"Added {added} new precursors. Curated file: {curated_faa}")
    if args.merge:
        print(f"Merged verified file: {out_dir / args.merge_output}")
        print(f"Merged verified TSV: {out_dir / args.merge_tsv_output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
