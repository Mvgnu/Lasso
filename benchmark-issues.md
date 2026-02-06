# Benchmark Issues (Pre-Presentation Note)

Date: 2026-02-05  
Scope: existing benchmark outputs are kept as-is for presentation; no reruns planned before presentation.

## Purpose
This note records known benchmark script/reporting caveats so interpretation is transparent.

## Known Issues

1. `loci_with_validated_precursor` naming is misleading in lasso holdout summary.
- Script: `scripts/benchmark_pipeline_validation.py`
- The value is incremented per evaluated locus-like event, not strict unique GBK files.
- Impact: label/readability issue.
- Performance bias risk: does **not** inflate top-k recall metrics.

2. Redundant summary fields in lasso holdout summary.
- Script: `scripts/benchmark_pipeline_validation.py`
- `topk` and `locus_level_topk` currently report the same calculation.
- Impact: presentation ambiguity (duplicate stats), not numerical corruption.
- Performance bias risk: none (same values, not optimistic values).

3. Silent exception handling while indexing GBK features.
- Scripts:
  - `scripts/benchmark_pipeline_validation.py`
  - `beta-lactamase-bench/scripts/benchmark_beta_lactamase.py`
- Some parse/record failures may be skipped without surfaced counts.
- Impact: potential undercount of evaluable loci.
- Performance bias risk: typically conservative/neutral; it does not create artificial “better” hits.

4. `record_id` alone in `holdout_topk.tsv` is not globally unique.
- Scripts:
  - `scripts/benchmark_pipeline_validation.py`
  - `beta-lactamase-bench/scripts/benchmark_beta_lactamase.py`
- Same `record_id` could appear across different GBK files.
- Impact: traceability issue when auditing rows manually.
- Performance bias risk: none for computed aggregate recall.

5. Core benchmark top-k helper is slightly unclear in denominator intent.
- Script: `scripts/benchmark_core_prediction_lab_dataset.py`
- Implementation treats unresolved ranks (`None`) as failures via denominator choice, but code style is easy to misread.
- Impact: readability/maintainability.
- Performance bias risk: conservative if anything (not optimistic).

## Interpretation for Presentation

- Reported benchmark recall values are still valid for ranking performance discussion.
- Known issues are primarily naming/reporting/traceability quality problems.
- None of the listed issues are known to artificially improve benchmark recall or create optimistic bias.

## Planned Post-Presentation Cleanup

- Align summary field names with true counting semantics.
- Remove redundant summary blocks.
- Log and count skipped records/exceptions explicitly.
- Add globally unique row identifiers (`gbk_path + record_id`) in top-k outputs.
- Clarify/standardize denominator semantics in core benchmark helper.
