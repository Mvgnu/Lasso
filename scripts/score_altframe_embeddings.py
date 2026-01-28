#!/usr/bin/env python3
"""
Score conserved alt-frame ORFs with existing ESM embedding + top-N mean cosine similarity.
License: MIT
Author: Magnus Ohle

Steps:
1) Filter conserved sequences from alt_frame_conservation summaries.
2) Reconstruct peptide sequences from GBKs using alt_frame_hits.tsv.
3) Embed candidates and score vs validated precursors.
4) Write altframe_conserved_candidates.faa + altframe_embedding_scores.tsv (+ nulls).
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import logging
import random
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd
import numpy as np
from Bio import SeqIO
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize

from lasso_workbench.pipeline.embed_precursor_candidates import (
    load_fasta,
    embed_records,
)
from lasso_workbench.pipeline.esm_embedder import ESM2Embedder
from lasso_workbench.pipeline.orf_extraction import chunk_orfs
from lasso_workbench.utils.sequence_io import write_fasta_pairs
from lasso_workbench.utils.translation import translate_bacterial

logger = logging.getLogger(__name__)


def _seq_hash(seq: str) -> str:
    return hashlib.sha1(seq.encode("utf-8")).hexdigest()[:12]


def _load_summary(
    summary_path: Path,
    min_genomes: int,
    require_conserved: bool,
    genes: Optional[set],
) -> pd.DataFrame:
    df = pd.read_csv(summary_path, sep="\t")
    if require_conserved and "conserved" in df.columns:
        conserved = df["conserved"].astype(str).str.lower().isin({"true", "1", "yes"})
        df = df[conserved]
    if min_genomes and "genomes" in df.columns:
        df = df[df["genomes"] >= min_genomes]
    if genes and "gene_name" in df.columns:
        df = df[df["gene_name"].isin(genes)]
    return df


def _sort_summary(df: pd.DataFrame) -> pd.DataFrame:
    sort_cols = [c for c in ("genomes", "hits") if c in df.columns]
    if sort_cols:
        return df.sort_values(by=sort_cols, ascending=False)
    return df


def _collect_hits_by_file(
    hits_path: Path,
    needed: set[str],
    max_per_hash: int = 1,
) -> Tuple[Dict[str, List[dict]], set[str], Dict[str, Dict[str, List[dict]]]]:
    per_file: Dict[str, List[dict]] = defaultdict(list)
    per_hash_counts: Dict[str, int] = defaultdict(int)
    found: set[str] = set()
    contexts: Dict[str, Dict[str, List[dict]]] = defaultdict(lambda: defaultdict(list))

    with hits_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            seq_hash = row.get("seq_hash")
            if not seq_hash or seq_hash not in needed:
                continue
            if per_hash_counts[seq_hash] >= max_per_hash:
                continue
            per_hash_counts[seq_hash] += 1
            found.add(seq_hash)
            gbk_file = row["gbk_file"]
            record_id = row["record_id"]
            per_file[gbk_file].append({
                "record_id": record_id,
                "orf_start": int(row["orf_start"]),
                "orf_end": int(row["orf_end"]),
                "orf_strand": row["orf_strand"],
                "seq_hash": seq_hash,
            })
            contexts[gbk_file][record_id].append({
                "seq_hash": seq_hash,
                "gene_start": int(row["gene_start"]),
                "gene_end": int(row["gene_end"]),
                "gene_strand": row.get("gene_strand", ""),
            })
            if len(found) == len(needed):
                break

    return per_file, found, contexts


def _extract_sequences(
    hits_by_file: Dict[str, List[dict]],
    gbk_dir: Path,
    needed: set[str],
) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    remaining = set(needed)

    for gbk_file, rows in hits_by_file.items():
        path = gbk_dir / gbk_file
        if not path.exists():
            logger.warning("Missing GBK file: %s", path)
            continue

        rows_by_record: Dict[str, List[dict]] = defaultdict(list)
        for row in rows:
            rows_by_record[row["record_id"]].append(row)

        try:
            for record_index, record in enumerate(SeqIO.parse(str(path), "genbank"), start=1):
                record_id = record.id or f"{path.stem}_record{record_index}"
                if record_id not in rows_by_record:
                    continue
                seq = record.seq
                for row in rows_by_record[record_id]:
                    seq_hash = row["seq_hash"]
                    if seq_hash in seqs:
                        continue
                    dna = seq[row["orf_start"]:row["orf_end"]]
                    if row["orf_strand"] == "-":
                        dna = dna.reverse_complement()
                    protein = translate_bacterial(str(dna), start_codon_to_met=True)
                    if not protein:
                        continue
                    computed = _seq_hash(protein)
                    if computed != seq_hash:
                        logger.warning("Hash mismatch for %s in %s", seq_hash, gbk_file)
                        continue
                    seqs[seq_hash] = protein
                    remaining.discard(seq_hash)
                if not remaining:
                    break
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed parsing %s: %s", path, exc)

        if not remaining:
            break

    if remaining:
        logger.warning("Missing %d sequences after GBK reconstruction", len(remaining))

    return seqs


def _build_candidates(seqs: Dict[str, str]) -> List[Tuple[str, str]]:
    pairs: List[Tuple[str, str]] = []
    for seq_hash, seq in sorted(seqs.items()):
        pairs.append((f"altframe|hash={seq_hash}", seq))
    return pairs


def _window_bounds(gene_start: int, gene_end: int, seq_len: int, pad_nt: int, max_aa: int) -> Tuple[int, int]:
    if pad_nt > 0:
        pad = pad_nt
    else:
        pad = max_aa * 3
    start = max(0, gene_start - pad)
    end = min(seq_len, gene_end + pad)
    if end <= start:
        return 0, seq_len
    return start, end


def _choose_random_orf(orfs, exclude_hash: str, rng: random.Random, max_tries: int = 12):
    if not orfs:
        return None
    for _ in range(max_tries):
        orf = rng.choice(orfs)
        if _seq_hash(orf.protein) != exclude_hash:
            return orf
    return None


def _sample_null_candidates(
    contexts: Dict[str, Dict[str, List[dict]]],
    sequences: Dict[str, str],
    gbk_dir: Path,
    min_aa: int,
    max_aa: int,
    gene_pad_nt: int,
    length_tol: int,
    per_candidate: int,
    seed: int,
) -> Tuple[List[Tuple[str, str]], pd.DataFrame]:
    rng = random.Random(seed)
    null_pairs: List[Tuple[str, str]] = []
    meta_rows: List[dict] = []
    missing = 0
    orf_cache: Dict[Tuple[str, int, int], list] = {}

    for gbk_file, record_map in contexts.items():
        path = gbk_dir / gbk_file
        if not path.exists():
            logger.warning("Missing GBK file for nulls: %s", path)
            continue
        try:
            for record_index, record in enumerate(SeqIO.parse(str(path), "genbank"), start=1):
                record_id = record.id or f"{path.stem}_record{record_index}"
                if record_id not in record_map:
                    continue
                for ctx in record_map[record_id]:
                    seq_hash = ctx["seq_hash"]
                    target_seq = sequences.get(seq_hash)
                    if not target_seq:
                        continue
                    target_len = len(target_seq)
                    window_start, window_end = _window_bounds(
                        ctx["gene_start"],
                        ctx["gene_end"],
                        len(record.seq),
                        gene_pad_nt,
                        max_aa,
                    )
                    cache_key = (record_id, window_start, window_end)
                    orfs = orf_cache.get(cache_key)
                    if orfs is None:
                        window_seq = record.seq[window_start:window_end]
                        forward = list(chunk_orfs(window_seq, "+", window_start, window_end, min_aa, max_aa))
                        rc_seq = window_seq.reverse_complement()
                        reverse = list(chunk_orfs(rc_seq, "-", window_start, window_end, min_aa, max_aa))
                        orfs = forward + reverse
                        orf_cache[cache_key] = orfs

                    if length_tol >= 0:
                        candidates = [
                            orf
                            for orf in orfs
                            if abs(int(orf.aa_len) - target_len) <= length_tol
                        ]
                    else:
                        candidates = [orf for orf in orfs if int(orf.aa_len) == target_len]
                    if not candidates:
                        missing += 1
                        continue
                    for sample_idx in range(per_candidate):
                        chosen = _choose_random_orf(candidates, seq_hash, rng)
                        if not chosen:
                            missing += 1
                            continue
                        protein = chosen.protein
                        null_hash = _seq_hash(protein)
                        candidate_id = f"null|source={seq_hash}|sample={sample_idx + 1}|hash={null_hash}"
                        null_pairs.append((candidate_id, protein))
                        meta_rows.append({
                            "candidate_id": candidate_id,
                            "source_seq_hash": seq_hash,
                            "null_seq_hash": null_hash,
                            "null_aa_length": len(protein),
                            "gbk_file": gbk_file,
                            "record_id": record_id,
                            "gene_start": ctx["gene_start"],
                            "gene_end": ctx["gene_end"],
                            "gene_strand": ctx.get("gene_strand", ""),
                        })
        except Exception as exc:  # noqa: BLE001
            logger.warning("Failed parsing %s for nulls: %s", path, exc)

    if missing:
        logger.warning("Null sampling missing for %d candidate windows", missing)

    return null_pairs, pd.DataFrame(meta_rows)


def _scores_to_df(rows: Sequence[dict]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    return pd.DataFrame(rows)


def _extract_locus_ids(reference_ids: Sequence[str]) -> List[str]:
    loci: List[str] = []
    for rid in reference_ids:
        locus = None
        for part in str(rid).split("|"):
            if part.startswith("locus="):
                locus = part.split("=", 1)[1]
                break
        if not locus:
            raise ValueError(f"Missing locus= token in reference id: {rid}")
        loci.append(locus)
    return loci


def _group_indices(loci: Sequence[str]) -> Tuple[List[str], List[List[int]]]:
    group_to_indices: Dict[str, List[int]] = {}
    for idx, locus in enumerate(loci):
        group_to_indices.setdefault(locus, []).append(idx)
    group_names = list(group_to_indices.keys())
    groups = list(group_to_indices.values())
    return group_names, groups


def _grouped_max(sim_matrix: np.ndarray, groups: List[List[int]]) -> np.ndarray:
    num_queries = sim_matrix.shape[0]
    grouped = np.empty((num_queries, len(groups)), dtype=np.float32)
    for group_idx, col_indices in enumerate(groups):
        grouped[:, group_idx] = sim_matrix[:, col_indices].max(axis=1)
    return grouped


def _topn_mean(grouped_matrix: np.ndarray, top_n: int) -> np.ndarray:
    num_groups = grouped_matrix.shape[1]
    top_n_count = max(1, int(top_n))
    if num_groups <= top_n_count:
        return grouped_matrix.mean(axis=1)
    partitioned = np.partition(grouped_matrix, -top_n_count, axis=1)
    top_vals = partitioned[:, -top_n_count:]
    return top_vals.mean(axis=1)


def _topk_indices_and_values(grouped_matrix: np.ndarray, k: int) -> Tuple[List[List[int]], List[List[float]]]:
    num_groups = grouped_matrix.shape[1]
    top_k = min(max(1, int(k)), num_groups)
    indices_list: List[List[int]] = []
    values_list: List[List[float]] = []
    for row in grouped_matrix:
        if num_groups <= top_k:
            idxs = np.argsort(row)[::-1]
        else:
            idxs = np.argpartition(row, -top_k)[-top_k:]
            idxs = idxs[np.argsort(row[idxs])[::-1]]
        indices_list.append(idxs.tolist())
        values_list.append(row[idxs].tolist())
    return indices_list, values_list


def _score_candidates_with_topk(
    validated_sequences: Sequence[Tuple[str, str]],
    candidate_sequences: Sequence[Tuple[str, str]],
    embedder: ESM2Embedder,
    score_batch_size: int,
    embed_batch_size: int,
    top_n_mean: int,
    top_k: int,
) -> List[dict]:
    ref_ids, ref_matrix = embed_records(embedder, validated_sequences, embed_batch_size=embed_batch_size)
    ref_norm = normalize(ref_matrix, norm="l2", axis=1, copy=True)
    locus_ids = _extract_locus_ids(ref_ids)
    group_names, groups = _group_indices(locus_ids)
    ref_ids_array = np.asarray(ref_ids)

    similarity_rows: List[dict] = []
    chunk_size = score_batch_size

    for start in range(0, len(candidate_sequences), chunk_size):
        chunk = candidate_sequences[start : start + chunk_size]
        cand_ids, cand_matrix = embed_records(embedder, chunk, embed_batch_size=embed_batch_size)
        cand_norm = normalize(cand_matrix, norm="l2", axis=1, copy=True)
        sim_matrix = cosine_similarity(cand_norm, ref_norm, dense_output=True)
        best_idx = sim_matrix.argmax(axis=1)
        best_scores = sim_matrix[np.arange(sim_matrix.shape[0]), best_idx]

        grouped_matrix = _grouped_max(sim_matrix, groups)
        top_n_means = _topn_mean(grouped_matrix, top_n_mean)
        topk_idxs, topk_vals = _topk_indices_and_values(grouped_matrix, top_k)

        for i, cand_id in enumerate(cand_ids):
            locus_ids = [group_names[idx] for idx in topk_idxs[i]]
            top_vals = topk_vals[i]
            similarity_rows.append({
                "candidate_id": cand_id,
                "best_match_id": str(ref_ids_array[best_idx[i]]),
                "best_similarity": float(best_scores[i]),
                "top_n_mean_similarity": float(top_n_means[i]),
                "top5_locus_sims": ",".join(f"{v:.4f}" for v in top_vals),
                "top5_locus_ids": ",".join(str(v) for v in locus_ids),
                "top5_min": float(min(top_vals)) if top_vals else 0.0,
            })

    return similarity_rows


def main() -> int:
    parser = argparse.ArgumentParser(description="Score conserved alt-frame ORFs with ESM embeddings.")
    parser.add_argument("--altframe-dir", type=Path, default=Path("results/alt_frame_conservation"))
    parser.add_argument("--summary", choices=["plain", "geom"], default="geom")
    parser.add_argument("--min-genomes", type=int, default=20)
    parser.add_argument("--max-seqs", type=int, default=500)
    parser.add_argument("--require-conserved", action="store_true", default=True)
    parser.add_argument("--no-require-conserved", action="store_false", dest="require_conserved")
    parser.add_argument("--genes", type=str, default="", help="Comma-separated gene names to keep")
    parser.add_argument("--hits-file", type=Path, default=None)
    parser.add_argument("--gbk-dir", type=Path, default=Path("data/antismash_lasso/gbk"))
    parser.add_argument("--validated-faa", type=Path, default=Path("data/precursors/precursor_proteins_verified.faa"))
    parser.add_argument("--model-name", type=str, default="facebook/esm2_t6_8M_UR50D")
    parser.add_argument("--device", choices=["auto", "cpu", "cuda", "mps"], default="auto")
    parser.add_argument("--min-aa", type=int, default=20)
    parser.add_argument("--max-aa", type=int, default=120)
    parser.add_argument("--embed-batch-size", type=int, default=100)
    parser.add_argument("--score-batch-size", type=int, default=100)
    parser.add_argument("--top-n-mean", type=int, default=5)
    parser.add_argument("--top-n-output", type=int, default=0, help="Keep only top N scored candidates (0 keeps all).")
    parser.add_argument(
        "--score-mode",
        choices=["top_n_mean", "best_similarity"],
        default="top_n_mean",
        help="Score column used for top-N output filtering.",
    )
    parser.add_argument("--gene-pad-nt", type=int, default=360)
    parser.add_argument("--null-per-candidate", type=int, default=1)
    parser.add_argument("--null-length-tol", type=int, default=0)
    parser.add_argument("--null-seed", type=int, default=13)
    parser.add_argument("--output-dir", type=Path, default=Path("results/altframe_embedding_scores"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    summary_path = args.altframe_dir / (
        "summary_sequences_geom.tsv" if args.summary == "geom" else "summary_sequences.tsv"
    )
    hits_path = args.hits_file or (args.altframe_dir / "alt_frame_hits.tsv")
    if not summary_path.exists():
        logger.error("Summary file not found: %s", summary_path)
        return 2
    if not hits_path.exists():
        logger.error("Hits file not found: %s", hits_path)
        return 2

    genes = {g.strip() for g in args.genes.split(",") if g.strip()} if args.genes else None
    summary_df = _load_summary(summary_path, args.min_genomes, args.require_conserved, genes)
    if summary_df.empty:
        logger.error("No conserved sequences after filtering.")
        return 2

    summary_df = _sort_summary(summary_df)
    if args.max_seqs and len(summary_df) > args.max_seqs:
        summary_df = summary_df.head(args.max_seqs)

    needed_hashes = set(summary_df["seq_hash"].dropna().astype(str).tolist())
    if not needed_hashes:
        logger.error("No seq_hash values in filtered summary.")
        return 2

    logger.info("Target seq_hash count: %d", len(needed_hashes))
    hits_by_file, found_hashes, contexts = _collect_hits_by_file(hits_path, needed_hashes, max_per_hash=1)
    if len(found_hashes) < len(needed_hashes):
        logger.warning("Only found %d/%d hashes in hits file", len(found_hashes), len(needed_hashes))

    sequences = _extract_sequences(hits_by_file, args.gbk_dir, needed_hashes)
    if not sequences:
        logger.error("No sequences reconstructed from GBKs.")
        return 2

    # Build candidates and write FASTA
    args.output_dir.mkdir(parents=True, exist_ok=True)
    candidates = _build_candidates(sequences)
    fasta_path = args.output_dir / "altframe_conserved_candidates.faa"
    write_fasta_pairs(candidates, fasta_path)
    logger.info("Wrote %d candidates to %s", len(candidates), fasta_path)

    null_pairs: List[Tuple[str, str]] = []
    null_meta = pd.DataFrame()
    if args.null_per_candidate > 0:
        null_pairs, null_meta = _sample_null_candidates(
            contexts=contexts,
            sequences=sequences,
            gbk_dir=args.gbk_dir,
            min_aa=args.min_aa,
            max_aa=args.max_aa,
            gene_pad_nt=args.gene_pad_nt,
            length_tol=args.null_length_tol,
            per_candidate=args.null_per_candidate,
            seed=args.null_seed,
        )
        if null_pairs:
            null_fasta = args.output_dir / "altframe_conserved_candidates_null.faa"
            write_fasta_pairs(null_pairs, null_fasta)
            logger.info("Wrote %d null candidates to %s", len(null_pairs), null_fasta)
        else:
            logger.warning("No null candidates generated.")

    validated = load_fasta(args.validated_faa)
    if not validated:
        logger.error("Validated FASTA empty or missing: %s", args.validated_faa)
        return 2

    device = None if args.device == "auto" else args.device
    embedder = ESM2Embedder(model_name=args.model_name, device=device)

    all_candidates = candidates + null_pairs
    similarity_rows = _score_candidates_with_topk(
        validated_sequences=validated,
        candidate_sequences=all_candidates,
        embedder=embedder,
        score_batch_size=args.score_batch_size,
        embed_batch_size=args.embed_batch_size,
        top_n_mean=args.top_n_mean,
        top_k=5,
    )

    scores_df = _scores_to_df(similarity_rows)
    if scores_df.empty:
        logger.error("No embedding scores produced.")
        return 2

    alt_mask = scores_df["candidate_id"].str.startswith("altframe|")
    null_mask = scores_df["candidate_id"].str.startswith("null|")

    seq_map = {f"altframe|hash={h}": h for h in sequences}
    scores_df.loc[alt_mask, "seq_hash"] = scores_df.loc[alt_mask, "candidate_id"].map(seq_map)
    summary_df = summary_df[summary_df["seq_hash"].isin(sequences.keys())].copy()
    summary_df["candidate_id"] = summary_df["seq_hash"].map(lambda h: f"altframe|hash={h}")
    summary_df["aa_length"] = summary_df["seq_hash"].map(lambda h: len(sequences[h]))
    summary_df["protein_sequence"] = summary_df["seq_hash"].map(sequences)

    output_df = summary_df.merge(scores_df[alt_mask], on="candidate_id", how="left")
    if args.top_n_output and args.top_n_output > 0:
        score_col = "top_n_mean_similarity" if args.score_mode == "top_n_mean" else "best_similarity"
        if score_col in output_df.columns:
            output_df = output_df.sort_values(score_col, ascending=False).head(args.top_n_output)
    out_path = args.output_dir / "altframe_embedding_scores.tsv"
    output_df.to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote embedding scores to %s", out_path)

    if null_pairs and not null_meta.empty:
        null_scores = scores_df[null_mask].merge(null_meta, on="candidate_id", how="left")
        null_path = args.output_dir / "altframe_embedding_scores_null.tsv"
        null_scores.to_csv(null_path, sep="\t", index=False)
        logger.info("Wrote null embedding scores to %s", null_path)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
