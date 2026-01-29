---
name: benchmarks
description: "Benchmark scripts for ranking/recall performance"
paths:
  scripts:
    - scripts/benchmark_pipeline_validation.py
    - scripts/benchmark_core_prediction_lab_dataset.py
    - beta-lactamase-bench/scripts/benchmark_beta_lactamase.py
    - beta-lactamase-bench/scripts/parse_ambler_table.py
    - beta-lactamase-bench/scripts/fetch_ambler_fastas.py
  data:
    - data/genomes
    - data/precursors/precursor_proteins_multi_strict.faa
    - data/precursors/lab_core_candidates_multi_strict.tsv
    - beta-lactamase-bench/data/ambler_class_a_beta_lactamases.faa
exports:
  - holdout_topk.tsv
  - holdout_candidates.tsv
  - holdout_summary.json
  - core_prediction_lab_summary.json
  - core_prediction_lab_cases.tsv
consumes:
  - pipeline.semantic_pipeline
  - pipeline.orf_extraction
  - core.embedding_scoring
  - core.prediction
verification:
  smoke: ".venv/bin/python scripts/benchmark_pipeline_validation.py --help"
---

# Benchmarks Domain

The **benchmarks** domain contains scripts used to quantify retrieval/ranking quality in a
"needle in haystack" setting.

## Core Ideas

- **Oracle window**: Use a known protein_id to select a local genomic window.
- **Enumerate all ORFs** in the window (no CDS filtering).
- **Holdout by accession**: remove all references matching the held-out accession.
- **Rank by embedding similarity** and report **top‑k recall**.

## Lasso Holdout (scripts/benchmark_pipeline_validation.py)

- Uses `pep=` accessions in validated FASTA headers.
- Window defaults to ±20 kb around the peptidase CDS.
- Outputs per‑locus recall + summary (top1/top5/top10/top50) and candidate count stats.
- Only recovered holdout candidates are saved to `holdout_candidates.tsv`.

## Beta-lactamase Holdout (beta-lactamase-bench/scripts/benchmark_beta_lactamase.py)

- Uses `prot=` accessions (and `locus=`) in validated FASTA headers.
- Window defaults to ±20 kb around the matching CDS in `data/genomes`.
- Same output format as the lasso benchmark.

## Data Requirements

- Validated FASTA **must** include `locus=` tokens in headers for grouped scoring.
  - Lasso: `pep=...|locus=...`
  - Beta‑lactamase: `prot=...|locus=...`

## Outputs

- `holdout_topk.tsv`: per‑locus top‑k success + candidate counts.
- `holdout_candidates.tsv`: only the recovered holdout ORFs and their ranks.
- `holdout_summary.json`: aggregate recall + candidate count stats.

## Core Prediction Benchmark (scripts/benchmark_core_prediction_lab_dataset.py)

- Evaluates rule-based core prediction on lab multi-candidate dataset.
- Uses `CorePredictor` + rules to rank cores within each precursor.
- Outputs:
  - `core_prediction_lab_summary.json`
  - `core_prediction_lab_cases.tsv`
