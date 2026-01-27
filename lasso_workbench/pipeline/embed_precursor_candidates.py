#!/usr/bin/env python3
"""Embed validated precursors and harvested ORFs with ESM-2 and score cosine similarity."""

from __future__ import annotations

from pathlib import Path
from typing import List, Optional, Sequence, Tuple, Callable

import numpy as np
from lasso_workbench.utils.sequence_io import read_fasta_pairs

from lasso_workbench.pipeline.esm_embedder import ESM2Embedder
from lasso_workbench.core.embedding_scoring import calculate_grouped_top_n_means
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize


def load_fasta(path: Path) -> List[Tuple[str, str]]:
    return read_fasta_pairs(path)


def embed_records(
    embedder: ESM2Embedder,
    records: Sequence[Tuple[str, str]],
    embed_batch_size: int,
) -> Tuple[List[str], np.ndarray]:
    ids = [name for name, _ in records]
    seqs = [seq for _, seq in records]
    matrix = embedder.embed_batch(seqs, batch_size=embed_batch_size)
    return ids, matrix


def score_candidates_in_memory(
    validated_sequences: Sequence[Tuple[str, str]],
    candidate_sequences: Sequence[Tuple[str, str]],
    embedder: ESM2Embedder,
    score_batch_size: int = 100,
    embed_batch_size: int = 100,
    top_n_mean: int = 5,
    progress_callback: Optional[Callable[[str], None]] = None,
) -> List[dict]:
    """
    Core function to embed and score candidates against validated references.
    
    Args:
        validated_sequences: List of (id, seq) for references
        candidate_sequences: List of (id, seq) for candidates
        embedder: Initialized ESM2Embedder
        score_batch_size: Number of candidates to process/score at once
        embed_batch_size: Batch size for ESM-2 inference
        top_n_mean: N parameter for Top-N Mean score
        progress_callback: Optional logging callback
        
    Returns:
        similarity_rows: List of scoring dicts
    """
    def _log(msg: str):
        if progress_callback:
            progress_callback(msg)

    # 1. Embed validated references -> ref_matrix [num_refs, hidden_dim]
    ref_ids, ref_matrix = embed_records(embedder, validated_sequences, embed_batch_size=embed_batch_size)
    
    # 2. L2-normalize refs (for cosine similarity via dot product)
    ref_norm = normalize(ref_matrix, norm="l2", axis=1, copy=True)

    # 3. Prepare candidate scoring loop (chunking, top-N aggregation, vectorized ref id handling)
    similarity_rows: List[dict] = []
    chunk_size = score_batch_size
    top_n_count = top_n_mean
    ref_ids_array = np.asarray(ref_ids)

    cand_count = 0
    total_candidates = len(candidate_sequences)
    
    _log(f"Scoring {total_candidates} candidates against {len(ref_ids)} references...")

    for start in range(0, total_candidates, chunk_size):
        # slice candidates into chunks
        chunk = candidate_sequences[start : start + chunk_size]
        
        # Embed this chunk -> cand_matrix [chunk_size, hidden_dim] aligned with cand_ids
        cand_ids, cand_matrix = embed_records(embedder, chunk, embed_batch_size=embed_batch_size)

        # L2-normalize candidate embeddings so dot products correspond to cosine similarity 
        cand_norm = normalize(cand_matrix, norm="l2", axis=1, copy=True)
        
        # Compute cosine similarity between each candidate and each reference
        # sim_matrix shape: [num_candidates_in_chunk, num_references]
        sim_matrix = cosine_similarity(cand_norm, ref_norm, dense_output=True)
        
        # Grouped Top-N Scoring (locus-group max -> top-N mean)
        top_n_means = calculate_grouped_top_n_means(sim_matrix, ref_ids, top_n=top_n_count)

        # Best single reference match per candidate (argmax over reference axis)
        best_idx = sim_matrix.argmax(axis=1)
        best_scores = sim_matrix[np.arange(sim_matrix.shape[0]), best_idx]
        
        for i, cand_id in enumerate(cand_ids):
            similarity_rows.append({
                "candidate_id": cand_id,
                "best_match_id": str(ref_ids_array[best_idx[i]]),
                "best_similarity": float(best_scores[i]),
                "top_n_mean_similarity": float(top_n_means[i]),
            })
        
        cand_count += len(cand_ids)
        if cand_count % (chunk_size * 5) == 0:
            _log(f"Processed {cand_count}/{total_candidates} candidates...")
             
    return similarity_rows
