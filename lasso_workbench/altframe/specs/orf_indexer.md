# Lasso ORF Indexer (Essence)

This module builds a local SQLite index of **20–120 aa ORFs** extracted from
**antiSMASH lasso‑annotated regions**. The index supports fast conservation
queries (e.g., “present in ≥ N genomes”) without rerunning ORF extraction.

## Inputs
- `gbk_dir`: directory containing `.gbk` files (antiSMASH outputs).
- Only regions whose `product/products` contain `"lasso"` are scanned.
- Only ORFs in the 20–120 aa range are indexed.

## Outputs (SQLite)
### Table: `peptides`
| column | type | description |
|---|---|---|
| `aa_hash` | TEXT (PK) | Blake2b hash of peptide sequence |
| `aa_seq` | TEXT | peptide sequence |
| `aa_len` | INTEGER | peptide length |

### Table: `peptide_counts`
| column | type | description |
|---|---|---|
| `aa_hash` | TEXT (PK) | FK to `peptides` |
| `instance_count` | INTEGER | total occurrences |

### Table: `peptide_genomes`
| column | type | description |
|---|---|---|
| `aa_hash` | TEXT | FK to `peptides` |
| `gbk_file` | TEXT | genome identifier |

### Table: `peptide_stats`
Materialized stats for fast UI queries:
| column | type | description |
|---|---|---|
| `aa_hash` | TEXT | FK to `peptides` |
| `aa_len` | INTEGER | peptide length |
| `genome_count` | INTEGER | distinct `gbk_file` count |
| `instance_count` | INTEGER | total instances |

## Usage (Python)
```python
from pathlib import Path
from lasso_workbench.altframe import build_lasso_orf_index

build_lasso_orf_index(
    gbk_dir=Path("data/antismash_lasso/gbk"),
    db_path=Path("results/lasso_orf_index.sqlite"),
    min_aa=20,
    max_aa=120,
    reset=True,
)
```

## Design Notes
- Uses the same ORF extraction logic as the main pipeline (`chunk_orfs`).
- Coordinates are 0‑based, end‑exclusive.
- Hashing uses Blake2b (16‑byte digest) for compact, collision‑resistant IDs.
