#!/usr/bin/env python3
import argparse
import csv
import json
import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional
from urllib.request import Request, urlopen
from urllib.error import HTTPError


API_URL = "https://antismash-db.secondarymetabolites.org/api/search"


def build_payload(
    term: str,
    offset: int,
    paginate: int,
    base_payload: Optional[Dict[str, object]] = None,
) -> Dict[str, object]:
    if base_payload:
        payload = json.loads(json.dumps(base_payload))
        payload["offset"] = offset
        payload["paginate"] = paginate
        try:
            payload["query"]["terms"]["value"] = term
        except Exception:
            payload["term"] = term
        return payload

    return {
        "query": {
            "search": "region",
            "terms": {
                "termType": "expr",
                "category": "type",
                "value": term,
                "filters": [],
                "count": 1,
            },
            "return_type": "json",
        },
        "paginate": paginate,
        "offset": offset,
    }


def fetch_page(payload: Dict[str, object]) -> Dict[str, object]:
    body = json.dumps(payload).encode("utf-8")
    req = Request(API_URL, data=body, method="POST")
    req.add_header("Content-Type", "application/json")
    req.add_header("Accept", "application/json")
    req.add_header("User-Agent", "lasso-workbench-antismash-scraper/1.0")
    try:
        with urlopen(req, timeout=60) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except HTTPError as exc:
        detail = ""
        try:
            detail = exc.read().decode("utf-8")
        except Exception:
            detail = ""
        raise RuntimeError(f"HTTP {exc.code}: {exc.reason}. Response: {detail}") from exc


def region_to_row(region: Dict[str, object]) -> Dict[str, object]:
    assembly_id = str(region.get("assembly_id", ""))
    record_number = region.get("record_number")
    region_number = region.get("region_number")
    cluster_url = ""
    gbk_url = ""
    if assembly_id and record_number is not None and region_number is not None:
        cluster_url = (
            f"https://antismash-db.secondarymetabolites.org/output/"
            f"{assembly_id}/index.html#r{record_number}c{region_number}"
        )
        gbk_url = (
            f"https://antismash-db.secondarymetabolites.org/output/"
            f"{assembly_id}/{assembly_id}.gbk"
        )
    mibig_acc = region.get("best_mibig_hit_acc") or ""
    mibig_url = f"https://mibig.secondarymetabolites.org/go/{mibig_acc}" if mibig_acc else ""

    return {
        "assembly_id": assembly_id,
        "acc": region.get("acc", ""),
        "record_number": record_number,
        "region_number": region_number,
        "start_pos": region.get("start_pos"),
        "end_pos": region.get("end_pos"),
        "description": region.get("description", ""),
        "cluster_url": cluster_url,
        "gbk_url": gbk_url,
        "best_mibig_hit_acc": mibig_acc,
        "best_mibig_hit_similarity": region.get("best_mibig_hit_similarity") or "",
        "mibig_url": mibig_url,
    }


def iter_regions(
    term: str,
    paginate: int,
    sleep_s: float,
    max_pages: int,
    base_payload: Optional[Dict[str, object]],
) -> Iterable[Dict[str, object]]:
    offset = 0
    pages = 0
    while True:
        payload = build_payload(term, offset, paginate, base_payload=base_payload)
        data = fetch_page(payload)
        regions: List[Dict[str, object]] = data.get("regions", [])
        if not regions:
            break
        for region in regions:
            yield region_to_row(region)
        offset += paginate
        pages += 1
        if max_pages > 0 and pages >= max_pages:
            break
        time.sleep(sleep_s)


def main() -> None:
    parser = argparse.ArgumentParser(description="Query antiSMASH DB API and export regions to CSV.")
    parser.add_argument("--term", default="lassopeptide", help="Search term")
    parser.add_argument("--paginate", type=int, default=50, help="Page size")
    parser.add_argument("--sleep", type=float, default=2.5, help="Seconds between requests")
    parser.add_argument("--max_pages", type=int, default=0, help="Limit pages (0 = all)")
    parser.add_argument(
        "--payload_json",
        default="",
        help="Optional JSON file containing the exact API payload (offset/paginate will be overwritten)",
    )
    parser.add_argument(
        "--payload_inline",
        default="",
        help="Optional JSON string for the exact API payload (offset/paginate will be overwritten)",
    )
    parser.add_argument("--out_csv", required=True, help="Output CSV path")
    args = parser.parse_args()

    base_payload = None
    if args.payload_json:
        base_payload = json.loads(Path(args.payload_json).read_text())
    elif args.payload_inline:
        base_payload = json.loads(args.payload_inline)

    rows = list(iter_regions(args.term, args.paginate, args.sleep, args.max_pages, base_payload))
    out_path = Path(args.out_csv)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    fieldnames = [
        "assembly_id",
        "acc",
        "record_number",
        "region_number",
        "start_pos",
        "end_pos",
        "description",
        "cluster_url",
        "gbk_url",
        "best_mibig_hit_acc",
        "best_mibig_hit_similarity",
        "mibig_url",
    ]
    with out_path.open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {out_path}")


if __name__ == "__main__":
    main()
