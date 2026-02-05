# Altframe Constraint Conservation

This package implements a constraint-aware conservation test for alt-frame ORFs
overlapping annotated CDS genes. It defines candidates by locus geometry
(gene + overlap bins + frame), then evaluates whether alt-frame conservation
exceeds what is expected from synonymous codon variation in the primary CDS.

## Inputs
- Hits file (`alt_frame_hits.tsv`):
  - Required columns:
    - `gene_name`, `match_type`, `orf_strand`, `orf_frame`
    - `gene_start`, `gene_end`, `orf_start`, `orf_end`
    - `gbk_file`, `record_id`
- GenBank directory:
  - GBK files with CDS features and `gene`/`locus_tag`/`product` annotations.

## Outputs
- `altframe_constraint_conservation.tsv` in the output directory with:
  - Locus keys (gene, match_type, strand, frame, bins)
  - Observed survival/identity
  - Null mean/std for survival/identity
  - z-score and p-values
  - Explicit count columns:
    - `hit_genomes`, `hit_events`: counts derived from the hits file
    - `evaluated_records`, `evaluated_instances`: counts used in observed/null scoring

## Core Concepts
### Coordinate convention
All coordinates are **0-based, end-exclusive** (`[start, end)`), matching
Biopython `FeatureLocation` and Python slicing. Bin assignment uses the last
covered base (`end - 1`) for `bin_end`.

### Frame convention
Canonical frame numbering matches the pipeline:
- Forward strand: `0/1/2`
- Reverse strand: `-1/-2/-3`

Use `--frame-mode one` only for legacy hits files that encode forward frames as
`1/2/3`. `--frame-mode auto` will infer this, but `zero` is the default.

### Locus definition
Each locus is keyed by:
`(gene_name, match_type, orf_strand, orf_frame, bin_start, bin_end)`.

The bins are computed from the overlap of an ORF with the gene and then
converted back into a window in gene-relative coordinates.

### Observed metrics
- **Survival**: fraction of genomes where the alt-frame window has no stop codon.
- **Identity**: mean pairwise identity among surviving peptides.

Implementation note:
- Scoring currently runs in **all-instances** mode for each selected locus key:
  all extracted gene instances for that gene are evaluated, not only the original
  hit rows. Use `evaluated_*` columns for denominator interpretation.

### Null model (synonymous shuffling)
For each genome, codons in the primary CDS are shuffled among synonymous
positions (codon 0 fixed) to preserve amino acid sequence and per-AA codon
multiset. The alt-frame peptide is recomputed from the same locus window.
If a CDS has `codon_start` != 1, the prefix bases are preserved and only the
in-frame codons are shuffled.

## Module Map
- `models.py`: `GeneInstance`, `WindowInstance`, `LocusKey`
- `binning.py`: bin computation + window extraction
- `conservation.py`: translation, identity, p-values
- `codon_shuffle.py`: synonymous shuffling + null sampling
- `gene_extraction.py`: GBK parsing and locus selection
- `cli.py`: orchestration + TSV output

## CLI
Entry point (script wrapper):
```
scripts/altframe_constraint_conservation.py
```

Example:
```
.venv/bin/python scripts/altframe_constraint_conservation.py \
  --altframe-dir results/alt_frame_conservation \
  --gbk-dir data/antismash_lasso/gbk \
  --min-genomes 20 \
  --max-candidates 200 \
  --geom-bins 20 \
  --null-iterations 200 \
  --output-dir results/altframe_constraint_tests
```

## Tests
Unit tests live in:
```
lasso_workbench/altframe/tests
```
They cover binning, translation/identity, shuffling, GBK extraction, and
CLI smoke execution.
