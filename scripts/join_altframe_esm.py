#!/usr/bin/env python3
"""
Join conserved alt-frame ORFs to ESM candidate/precursor sets (exact or k-mer).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import logging
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq

from lasso_workbench.utils.translation import translate_bacterial

logger = logging.getLogger(__name__)


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:12]


def _load_summary(summary_path: Path, min_genomes: int, require_conserved: bool, genes: Optional[set]) -> pd.DataFrame:
    df = pd.read_csv(summary_path, sep="	")
    if require_conserved and "conserved" in df.columns:
        df = df[df["conserved"] == True]  # noqa: E712
    if min_genomes and "genomes" in df.columns:
        df = df[df["genomes"] >= min_genomes]
    if genes and "gene_name" in df.columns:
        df = df[df["gene_name"].isin(genes)]
    return df


def _collect_needed_hashes(summary_df: pd.DataFrame) -> set[str]:
    return set(summary_df["seq_hash"].dropna().astype(str).tolist())


def _group_hits_by_file(hits_path: Path, needed: set[str]) -> Dict[str, List[dict]]:
    per_file: Dict[str, List[dict]] = defaultdict(list)
    with hits_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="	")
        for row in reader:
            seq_hash = row.get("seq_hash")
            if not seq_hash or seq_hash not in needed:
                continue
            per_file[row["gbk_file"]].append({
                "record_id": row["record_id"],
                "orf_start": int(row["orf_start"]),
                "orf_end": int(row["orf_end"]),
                "orf_strand": row["orf_strand"],
                "seq_hash": seq_hash,
            })
    return per_file


def _extract_sequences(
    hits_by_file: Dict[str, List[dict]],
    gbk_dir: Path,
    needed: set[str],
) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    missing_files = 0
    for gbk_file, rows in hits_by_file.items():
        path = gbk_dir / gbk_file
        if not path.exists():
            missing_files += 1
            logger.warning("Missing GBK file: %s", path)
            continue
        rows_by_record: Dict[str, List[dict]] = defaultdict(list)
        for row in rows:
            rows_by_record[row["record_id"]].append(row)
        for record_index, record in enumerate(SeqIO.parse(str(path), "genbank"), start=1):
            record_id = record.id or f"{path.stem}_record{record_index}"
            if record_id not in rows_by_record:
                continue
            seq = record.seq
            for row in rows_by_record[record_id]:
                if row["seq_hash"] in seqs:
                    continue
                dna = seq[row["orf_start"]:row["orf_end"]]
                if row["orf_strand"] == "-":
                    dna = dna.reverse_complement()
                protein = translate_bacterial(dna, start_codon_to_met=True)
                if not protein:
                    continue
                computed = _seq_hash(protein)
                if computed != row["seq_hash"]:
                    logger.warning("Hash mismatch for %s in %s", row["seq_hash"], gbk_file)
                if row["seq_hash"] in needed and row["seq_hash"] not in seqs:
                    seqs[row["seq_hash"]] = protein
        if len(seqs) == len(needed):
            break
    if missing_files:
        logger.warning("Missing %d GBK files during sequence reconstruction", missing_files)
    return seqs


def _load_candidates_tsv(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="	")
    if "protein_sequence" not in df.columns:
        raise ValueError("protein_sequence column missing in candidates TSV")
    return df


def _load_fasta(path: Path) -> pd.DataFrame:
    rows = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq).strip().upper()
        if not seq:
            continue
        rows.append({"record_id": record.id, "protein_sequence": seq})
    return pd.DataFrame(rows)


def _exact_join(alt_df: pd.DataFrame, target_df: pd.DataFrame) -> pd.DataFrame:
    target = target_df.copy()
    target["seq_hash"] = target["protein_sequence"].astype(str).map(_seq_hash)
    joined = target.merge(alt_df, on="seq_hash", how="inner", suffixes=("_target", "_alt"))
    return joined


def _kmer_index(seqs: Dict[str, str], k: int) -> Dict[str, List[str]]:
    index: Dict[str, List[str]] = defaultdict(list)
    for seq_hash, seq in seqs.items():
        seq = seq.upper()
        if len(seq) < k:
            continue
        seen = set()
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer in seen:
                continue
            seen.add(kmer)
            index[kmer].append(seq_hash)
    return index


def _kmer_join(alt_df: pd.DataFrame, alt_seqs: Dict[str, str], target_df: pd.DataFrame, k: int, min_hits: int) -> pd.DataFrame:
    index = _kmer_index(alt_seqs, k)
    rows = []
    for _, row in target_df.iterrows():
        seq = str(row.get("protein_sequence", "")).upper()
        if len(seq) < k:
            continue
        counts: Dict[str, int] = defaultdict(int)
        seen = set()
        for i in range(0, len(seq) - k + 1):
            kmer = seq[i:i+k]
            if kmer in seen:
                continue
            seen.add(kmer)
            for seq_hash in index.get(kmer, []):
                counts[seq_hash] += 1
        for seq_hash, count in counts.items():
            if count >= min_hits:
                rows.append({
                    **row.to_dict(),
                    "seq_hash": seq_hash,
                    "kmer_hits": count,
                })
    if not rows:
        return pd.DataFrame()
    joined = pd.DataFrame(rows).merge(alt_df, on="seq_hash", how="inner", suffixes=("_target", "_alt"))
    return joined


def main() -> int:
    parser = argparse.ArgumentParser(description="Join conserved alt-frame ORFs to ESM candidates.")
    parser.add_argument("--altframe-dir", type=Path, default=Path("results/alt_frame_conservation"))
    parser.add_argument("--summary", choices=["plain", "geom"], default="geom")
    parser.add_argument("--min-genomes", type=int, default=10)
    parser.add_argument("--require-conserved", action="store_true")
    parser.add_argument("--genes", type=str, default="", help="Comma-separated gene names to keep")
    parser.add_argument("--hits-file", type=Path, default=None)
    parser.add_argument("--gbk-dir", type=Path, default=Path("data/antismash_lasso/gbk"))
    parser.add_argument("--candidates-tsv", type=Path, default=None)
    parser.add_argument("--precursors-fasta", type=Path, default=None)
    parser.add_argument("--match", choices=["exact", "kmer"], default="exact")
    parser.add_argument("--kmer", type=int, default=7)
    parser.add_argument("--min-kmer-hits", type=int, default=3)
    parser.add_argument("--output-dir", type=Path, default=Path("results/altframe_esm_join"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    summary_path = args.altframe_dir / ("summary_sequences_geom.tsv" if args.summary == "geom" else "summary_sequences.tsv")
    if not summary_path.exists():
        logger.error("Summary file not found: %s", summary_path)
        return 2

    genes = {g.strip() for g in args.genes.split(",") if g.strip()} if args.genes else None
    summary_df = _load_summary(summary_path, args.min_genomes, args.require_conserved, genes)
    if summary_df.empty:
        logger.error("No alt-frame sequences after filtering.")
        return 2

    needed = _collect_needed_hashes(summary_df)
    hits_path = args.hits_file or (args.altframe_dir / "alt_frame_hits.tsv")
    if not hits_path.exists():
        logger.error("Hits file not found: %s", hits_path)
        return 2

    logger.info("Filtering %d alt-frame seq_hashes", len(needed))
    hits_by_file = _group_hits_by_file(hits_path, needed)
    alt_seqs = _extract_sequences(hits_by_file, args.gbk_dir, needed)
    if not alt_seqs:
        logger.error("No sequences reconstructed from hits.")
        return 2

    alt_df = summary_df.copy()
    alt_df["protein_sequence"] = alt_df["seq_hash"].map(alt_seqs)
    alt_df = alt_df.dropna(subset=["protein_sequence"])

    args.output_dir.mkdir(parents=True, exist_ok=True)
    alt_out = args.output_dir / "altframe_sequences.tsv"
    alt_df.to_csv(alt_out, sep="	", index=False)
    logger.info("Wrote %s", alt_out)

    if args.candidates_tsv:
        cand_df = _load_candidates_tsv(args.candidates_tsv)
        if args.match == "exact":
            joined = _exact_join(alt_df, cand_df)
        else:
            joined = _kmer_join(alt_df, alt_seqs, cand_df, args.kmer, args.min_kmer_hits)
        out_path = args.output_dir / ("join_candidates_kmer.tsv" if args.match == "kmer" else "join_candidates_exact.tsv")
        joined.to_csv(out_path, sep="	", index=False)
        logger.info("Wrote %s", out_path)

    if args.precursors_fasta:
        prec_df = _load_fasta(args.precursors_fasta)
        if args.match == "exact":
            joined = _exact_join(alt_df, prec_df)
        else:
            joined = _kmer_join(alt_df, alt_seqs, prec_df, args.kmer, args.min_kmer_hits)
        out_path = args.output_dir / ("join_precursors_kmer.tsv" if args.match == "kmer" else "join_precursors_exact.tsv")
        joined.to_csv(out_path, sep="	", index=False)
        logger.info("Wrote %s", out_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
