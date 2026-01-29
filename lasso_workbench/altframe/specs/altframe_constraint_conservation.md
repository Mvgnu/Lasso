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

## Core Concepts
### Locus definition
Each locus is keyed by:
`(gene_name, match_type, orf_strand, orf_frame, bin_start, bin_end)`.

The bins are computed from the overlap of an ORF with the gene and then
converted back into a window in gene-relative coordinates.

### Observed metrics
- **Survival**: fraction of genomes where the alt-frame window has no stop codon.
- **Identity**: mean pairwise identity among surviving peptides.

### Null model (synonymous shuffling)
For each genome, codons in the primary CDS are shuffled among synonymous
positions (codon 0 fixed) to preserve amino acid sequence and per-AA codon
multiset. The alt-frame peptide is recomputed from the same locus window.

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
