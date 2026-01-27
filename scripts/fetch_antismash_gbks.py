#!/usr/bin/env python3
import argparse
import csv
import random
import time
from pathlib import Path
from urllib.request import urlretrieve


def main() -> None:
    parser = argparse.ArgumentParser(description="Download antiSMASH assembly GBKs from CSV.")
    parser.add_argument("--csv", required=True, help="CSV from antismash_api_search.py")
    parser.add_argument("--out_dir", required=True, help="Output directory for GBKs")
    parser.add_argument("--min_sleep", type=float, default=3.0, help="Minimum sleep between requests")
    parser.add_argument("--max_sleep", type=float, default=5.0, help="Maximum sleep between requests")
    parser.add_argument("--max_rows", type=int, default=0, help="Limit rows (0 = all)")
    parser.add_argument("--overwrite", action="store_true", help="Overwrite existing files")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    total = 0
    with open(args.csv, newline="") as f:
        total = sum(1 for _ in f) - 1
    if total < 0:
        total = 0
    if args.max_rows > 0:
        total = min(total, args.max_rows)

    processed = 0
    with open(args.csv, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            if args.max_rows > 0 and processed >= args.max_rows:
                break
            processed += 1
            idx = processed
            url = (row.get("gbk_url") or "").strip()
            assembly_id = (row.get("assembly_id") or "").strip()
            if not url or not assembly_id:
                print(f"[{idx}/{total}] Missing url/assembly_id, skipping")
                continue

            out_path = out_dir / f"{assembly_id}.gbk"
            if out_path.exists() and not args.overwrite:
                print(f"[{idx}/{total}] Exists, skipping {out_path.name}")
                continue

            try:
                print(f"[{idx}/{total}] Downloading {assembly_id}")
                urlretrieve(url, out_path)
            except Exception as e:
                print(f"[{idx}/{total}] Failed {url}: {e}")

            sleep_s = random.uniform(args.min_sleep, args.max_sleep)
            time.sleep(sleep_s)


if __name__ == "__main__":
    main()
