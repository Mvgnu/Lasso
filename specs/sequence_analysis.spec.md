---
name: sequence_analysis
description: "Core logic for sequence feature extraction, machinery detection, and analysis"
paths:
  modules:
    - lasso_workbench/core/embedding_scoring.py
    - lasso_workbench/pipeline/orf_extraction.py
    - lasso_workbench/pipeline/semantic_pipeline.py
  tests:
    - tests/test_scientific_features.py
    - tests/test_consolidation.py
exports:
  - calculate_similarity_matrix
  - extract_all_orfs
consumes:
  - bio.seq
  - numpy
  - pandas
  - sklearn.metrics.pairwise
verification:
  test: "pytest tests/test_scientific_features.py tests/test_consolidation.py -v"
---

# Sequence Analysis Domain

The **sequence_analysis** domain handles the extraction of biological features from raw sequence data. This includes ORF extraction and vector embedding analysis.

## Core Philosophy: "No Slop"

**Logic must be consolidated into specific libraries.** Do not write ad-hoc loops or manual string parsers in control flow.

> [!IMPORTANT]
> **Use Established Libraries**:
> *   **Sequences**: Use `BioPython` (`Bio.Seq`, `Bio.SeqIO`). Do not manually translate or parse FASTA.
> *   **Data**: Use `Pandas` for tabular data (candidates, scores). Avoid lists of dictionaries where DataFrames facilitate vectorized ops.
> *   **Math**: Use `NumPy` and `Scikit-learn` for vector operations (cosine similarity, normalization). Do not write manual dot products.

## Module Responsibilities

| Module | Responsibility | Library Enforced |
|--------|----------------|------------------|
| `core/embedding_scoring.py` | Vector scoring & normalization | `numpy`, `sklearn` |
| `pipeline/orf_extraction.py` | 6-frame translation & ORF filtering | `Bio.Seq` |
| `pipeline/semantic_pipeline.py` | Orchestration of analysis steps | `pandas` (merging results) |

## Anti-Patterns to Avoid

*   **Loop nesting for scoring**: Do not loop `for candidate in candidates: for ref in refs: score()`. Use matrix operations.
*   **Manual Translation**: Do not use dictionaries for codon tables. Use `Bio.Seq.translate(table=11)`.
*   **Inline Heuristics**: Do not define `def _check_motif(s):` inside a pipeline function. Move it to a core library (`core/prediction.py` or similar).

## Data Flow

1.  **Input**: Raw DNA (GBK/FASTA) → `Bio.SeqRecord`
2.  **Extraction**: `pipeline/orf_extraction` → `pd.DataFrame` (Candidates)
3.  **Enrichment**: `core/machinery` + `core/prediction` → Added Columns
4.  **Embedding**: `pipeline/embed_precursor_candidates` → `np.ndarray`
5.  **Scoring**: `core/embedding_scoring` → `pd.DataFrame` (Scores)

Keep the pipeline simple. Rely on the libraries.
