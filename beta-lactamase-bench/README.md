# Ambler Class A Beta-Lactamase Benchmark

Goal: run the same "needle in haystack" ranking benchmark on class A beta-lactamases.

## 1) Parse the Ambler table HTML

```bash
.venv/bin/python beta-lactamase-bench/scripts/parse_ambler_table.py \
  --table-html data/ambler-class-a-beta-lactamases/table.html \
  --output-tsv beta-lactamase-bench/data/ambler_class_a_table.tsv
```

## 2) Fetch FASTA sequences from NCBI

```bash
.venv/bin/python beta-lactamase-bench/scripts/fetch_ambler_fastas.py \
  --table-tsv beta-lactamase-bench/data/ambler_class_a_table.tsv \
  --output-faa beta-lactamase-bench/data/ambler_class_a_beta_lactamases.faa \
  --output-tsv beta-lactamase-bench/data/ambler_class_a_beta_lactamases.tsv \
  --email you@example.com
```

This creates headers like:
```
>prot=AAC09015|locus=AAC09015|name=AER-1|class=A|source=ambler_table
```

## What your data needs to look like

### FASTA headers (required)

The benchmark expects a **protein accession** and a **locus token** so embedding scores can be grouped:

```
>prot=ACCESSION|locus=ACCESSION|name=NAME|class=A|source=ambler_table
```

Notes:
- `prot=` is used for holdout selection.
- `locus=` is required by the embedding scorer (grouped top‑N mean).
- Any extra fields are fine.

### Genomes

Your genome folder must contain GBK/GenBank files with `CDS` features that include a matching `protein_id`.
Those `protein_id`s are used as the oracle to define the benchmark window.

## 3) Run the benchmark

```bash
.venv/bin/python beta-lactamase-bench/scripts/benchmark_beta_lactamase.py \
  --gbk_dir data/genomes \
  --validated_faa beta-lactamase-bench/data/ambler_class_a_beta_lactamases.faa \
  --out_dir beta-lactamase-bench/results/holdout_8m \
  --model_name facebook/esm2_t6_8M_UR50D \
  --device cpu
```

Recommended defaults:
- `--min_aa 200 --max_aa 400` (class A beta-lactamases are ~280 aa)
- `--window_nt 20000` (40 kb total window)

Outputs:
- `holdout_topk.tsv`
- `holdout_candidates.tsv` (only recovered holdouts)
- `holdout_summary.json`

Example summary (8M, 5 loci):
```
{
  "topk": {"top1": 0.8, "top5": 1.0, "top10": 1.0, "top50": 1.0},
  "candidate_count": {"mean": 354.6, "median": 395.0, "min": 276, "max": 395}
}
```
