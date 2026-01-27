"""
Core embedding scoring algorithms for Lasso peptide precursor discovery.
License: MIT
Author: Magnus Ohle
"""

from __future__ import annotations

from typing import List, Dict

import numpy as np
def calculate_grouped_top_n_means(
    similarity_matrix: np.ndarray,
    reference_ids: List[str],
    top_n: int = 5
) -> np.ndarray:
    """
     Compute a "Grouped Top-N Mean" score per query.

    References are grouped by their locus identifier extracted from the "locus=" identifier
    in the reference id. For each query and each locus-group, we take the maximum
    similarity across all reference variants that belong to that locus. This yields a
    per-query vector of locus-level best scores and avoids multiple sequences from the same locus 
    dominating the score.

    The final score per query is the mean of the top-N locus-level scores
    (or the mean across all loci if there are fewer than N loci).

    Args:
        similarity_matrix: [Q, R] cosine similarity matrix (Q queries x R references).
        reference_ids: Length-R list of reference identifiers (must contain "locus=" token).
        top_n: Number of top locus scores to average (clamped to >= 1).

    Returns:
        scores: [Q] array of grouped Top-N mean scores (one per query).
    """
    num_queries = similarity_matrix.shape[0]

    # Extract locus key per reference (requires locus=...)
    group_ids: List[str] = []
    for rid in reference_ids:
        locus = None
        for part in str(rid).split("|"):
            if part.startswith("locus="):
                locus = part.split("=", 1)[1]
                break
        if not locus:
            raise ValueError(f"Missing locus= token in reference id: {rid}")
        group_ids.append(locus)
    # 2) Collect column indices for each group key.
    group_to_indices: Dict[str, List[int]] = {}
    for idx, gid in enumerate(group_ids):
        group_to_indices.setdefault(gid, []).append(idx)

    # 3) Collapse reference variants within each locus by keeping only the per-query max
    groups = list(group_to_indices.values())
    num_groups = len(groups)
    grouped_matrix = np.empty((num_queries, num_groups), dtype=np.float32)
    for group_idx, col_indices in enumerate(groups):
        # shape: [num_queries], one max score per query for this locus/group
        grouped_matrix[:, group_idx] = similarity_matrix[:, col_indices].max(axis=1)
             
    # 4) Per query, average the top-N locus-max scores
    top_n_count = max(1, int(top_n))
    
    if num_groups <= top_n_count:
        # If fewer groups than N, take mean of all groups.
        return grouped_matrix.mean(axis=1)
    
    # Otherwise, use partition to find top N elements efficiently.
    # np.partition moves the k-th smallest element to position k.
    # We want top N, so we partition by -top_n_count.
    partitioned = np.partition(grouped_matrix, -top_n_count, axis=1)
    top_n_vals = partitioned[:, -top_n_count:]
    
    return top_n_vals.mean(axis=1)
