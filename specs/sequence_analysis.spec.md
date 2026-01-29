---
name: sequence_analysis
description: "Core logic for sequence feature extraction, machinery detection, and analysis"
paths:
  modules:
    - lasso_workbench/core/embedding_scoring.py
    - lasso_workbench/core/ranking.py
    - lasso_workbench/core/prediction.py
    - lasso_workbench/pipeline/orf_extraction.py
    - lasso_workbench/pipeline/embed_precursor_candidates.py
    - lasso_workbench/pipeline/esm_embedder.py
    - lasso_workbench/pipeline/semantic_pipeline.py
  tests:
    - tests/test_scientific_features.py
    - tests/test_consolidation.py
exports:
  - calculate_grouped_top_n_means
  - chunk_orfs
  - run_semantic_pipeline
  - score_candidates_in_memory
consumes:
  - bio.seq
  - numpy
  - pandas
  - sklearn.metrics.pairwise
  - torch
verification:
  test: "pytest tests/test_scientific_features.py tests/test_consolidation.py -v"
---

# Sequence Analysis Domain

The **sequence_analysis** domain handles the extraction of biological features from raw sequence data. This includes ORF extraction and vector embedding analysis.

## Core Philosophy: "No Slop"

**Logic must be consolidated into specific libraries.** Do not write ad-hoc loops or manual string parsers in control flow.

> [!IMPORTANT]
> **Use Established Libraries** (core pipeline only):
> *   **Sequences**: Use `BioPython` (`Bio.Seq`, `Bio.SeqIO`) for translation and GenBank parsing.
> *   **FASTA/TSV in scripts**: small scripts may parse FASTA/TSV manually for speed or to avoid heavy deps, but core modules should not.
> *   **Data**: Use `Pandas` for tabular data (candidates, scores). Avoid lists of dictionaries where DataFrames facilitate vectorized ops.
> *   **Math**: Use `NumPy` and `Scikit-learn` for vector operations (cosine similarity, normalization). Do not write manual dot products.

## Module Responsibilities

| Module | Responsibility | Library Enforced |
|--------|----------------|------------------|
| `core/embedding_scoring.py` | Vector scoring & normalization | `numpy`, `sklearn` |
| `core/ranking.py` | Embedding-first ranking + rule cutoff | - |
| `core/prediction.py` | Rule-engine core prediction | `BioPython` |
| `pipeline/orf_extraction.py` | 6-frame translation & ORF filtering | `Bio.Seq` |
| `pipeline/embed_precursor_candidates.py` | Embedding + cosine scoring | `numpy`, `sklearn` |
| `pipeline/esm_embedder.py` | ESM-2 embedding | `torch` |
| `pipeline/semantic_pipeline.py` | Orchestration of analysis steps | `pandas` (merging results) |

## Anti-Patterns to Avoid

*   **Loop nesting for scoring**: Do not loop `for candidate in candidates: for ref in refs: score()`. Use matrix operations.
*   **Manual Translation in core**: Do not use dictionaries for codon tables. Use `Bio.Seq.translate(table=11)` or `translate_bacterial`.
*   **Inline Heuristics**: Do not define `def _check_motif(s):` inside a pipeline function. Move it to a core library (`core/prediction.py` or similar).

## Data Flow

1.  **Input**: Raw DNA (GBK/FASTA) → `Bio.SeqRecord`
2.  **Extraction**: `pipeline/orf_extraction` → ORF candidates (6-frame)
3.  **Embedding**: `pipeline/embed_precursor_candidates` → `np.ndarray`
4.  **Scoring**: `core/embedding_scoring` → similarity rows
5.  **Ranking**: `core/ranking` applies optional rule cutoff (no weighted combo score)
6.  **Results**: `pipeline/semantic_pipeline` builds `PipelineResult` objects + TSV/JSON

Keep the pipeline simple. Rely on the libraries.
