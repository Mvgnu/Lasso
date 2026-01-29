#!/usr/bin/env python3
"""
Parse Ambler class A beta-lactamase HTML table into a TSV metadata file.
License: MIT
Author: Magnus Ohle
"""
from __future__ import annotations

import argparse
import csv
import html
from dataclasses import dataclass
from html.parser import HTMLParser
from pathlib import Path
from typing import List, Optional


@dataclass
class Cell:
    text: str
    href: str = ""


class AmblerTableParser(HTMLParser):
    def __init__(self) -> None:
        super().__init__()
        self.rows: List[List[Cell]] = []
        self._row: List[Cell] = []
        self._in_td = False
        self._cell_text: List[str] = []
        self._cell_href: str = ""

    def handle_starttag(self, tag: str, attrs) -> None:
        if tag == "tr":
            self._row = []
        elif tag == "td":
            self._in_td = True
            self._cell_text = []
            self._cell_href = ""
        elif tag == "a" and self._in_td:
            for key, value in attrs:
                if key == "href" and value:
                    if not self._cell_href:
                        self._cell_href = value
                    break

    def handle_endtag(self, tag: str) -> None:
        if tag == "td":
            text = html.unescape(" ".join(self._cell_text)).strip()
            self._row.append(Cell(text=text, href=self._cell_href))
            self._in_td = False
        elif tag == "tr":
            if self._row:
                self.rows.append(self._row)
            self._row = []

    def handle_data(self, data: str) -> None:
        if self._in_td:
            text = data.strip()
            if text:
                self._cell_text.append(text)


def normalize_accession(value: str) -> str:
    raw = str(value or "").strip()
    if not raw:
        return ""
    raw = raw.split()[0].strip()
    return raw.replace("\u00a0", "").strip()


def parse_table(html_path: Path) -> List[dict]:
    parser = AmblerTableParser()
    parser.feed(html_path.read_text())

    rows: List[dict] = []
    for cells in parser.rows:
        if len(cells) < 8:
            continue
        ambler_class = cells[0].text.strip()
        if ambler_class != "A":
            continue
        if len(cells) < 13:
            # Pad missing cells to keep column mapping stable
            cells = cells + [Cell("")] * (13 - len(cells))
        rows.append(
            {
                "ambler_class": ambler_class,
                "protein_name": cells[1].text.strip(),
                "alt_names": cells[2].text.strip(),
                "subfamily": cells[3].text.strip(),
                "genpept_id": normalize_accession(cells[4].text),
                "genpept_link": cells[4].href,
                "genbank_id": normalize_accession(cells[5].text),
                "genbank_link": cells[5].href,
                "pubmed": cells[6].text.strip(),
                "sequence_link": cells[7].href,
                "pdb_count": cells[8].text.strip(),
                "mutants": cells[9].text.strip(),
                "phenotype": cells[10].text.strip(),
                "functional_info": cells[11].text.strip(),
                "natural_or_acquired": cells[12].text.strip(),
            }
        )
    return rows


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--table-html",
        type=Path,
        default=Path("data/ambler-class-a-beta-lactamases/table.html"),
    )
    ap.add_argument(
        "--output-tsv",
        type=Path,
        default=Path("beta-lactamase-bench/data/ambler_class_a_table.tsv"),
    )
    args = ap.parse_args()

    rows = parse_table(args.table_html)
    if not rows:
        raise SystemExit("No class A entries found in table.")

    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {args.output_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
