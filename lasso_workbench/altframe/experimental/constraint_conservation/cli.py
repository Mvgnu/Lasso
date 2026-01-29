from __future__ import annotations

import argparse
import logging
import math
import random
from pathlib import Path
from typing import List

import pandas as pd

from lasso_workbench.altframe.experimental.constraint_conservation.binning import prepare_locus_windows
from lasso_workbench.altframe.experimental.constraint_conservation.codon_shuffle import null_for_locus
from lasso_workbench.altframe.experimental.constraint_conservation.conservation import observed_for_locus, p_value, summarize_null
from lasso_workbench.altframe.experimental.constraint_conservation.gene_extraction import (
    extract_gene_instances,
    scan_hits_for_loci,
    select_top_loci,
)

logger = logging.getLogger(__name__)


def main() -> int:
    parser = argparse.ArgumentParser(description="Constraint-aware conservation test for alt-frame loci.")
    parser.add_argument("--altframe-dir", type=Path, default=Path("results/alt_frame_conservation"))
    parser.add_argument("--min-genomes", type=int, default=20)
    parser.add_argument("--max-candidates", type=int, default=200)
    parser.add_argument("--geom-bins", type=int, default=20)
    parser.add_argument("--null-iterations", type=int, default=200)
    parser.add_argument("--seed", type=int, default=13)
    parser.add_argument("--gbk-dir", type=Path, default=Path("data/antismash_lasso/gbk"))
    parser.add_argument("--hits-file", type=Path, default=None)
    parser.add_argument(
        "--loci-file",
        type=Path,
        default=None,
        help="Optional TSV of loci to re-score (e.g., a previous output file).",
    )
    parser.add_argument(
        "--loci-limit",
        type=int,
        default=None,
        help="Optional limit on number of loci loaded from --loci-file.",
    )
    parser.add_argument("--gene-field", choices=["gene", "locus_tag", "product", "any"], default="gene")
    parser.add_argument(
        "--frame-mode",
        choices=["auto", "zero", "one"],
        default="zero",
        help=(
            "Frame numbering convention. Canonical pipeline frames are "
            "zero=+0/+1/+2 and -1/-2/-3. Use one=+1/+2/+3 only for legacy hits. "
            "auto infers from hits."
        ),
    )
    parser.add_argument("--output-dir", type=Path, default=Path("results/altframe_constraint_tests"))
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")

    hits_path = args.hits_file or (args.altframe_dir / "alt_frame_hits.tsv")

    loci: List[LocusKey] = []
    frame_values: set[int] = set()
    unique_counts = {}
    hit_counts = {}

    if args.loci_file:
        if not args.loci_file.exists():
            logger.error("Loci file not found: %s", args.loci_file)
            return 2
        loci_df = pd.read_csv(args.loci_file, sep="\t")
        required = {"gene_name", "match_type", "orf_strand", "orf_frame", "bin_start", "bin_end"}
        missing = required - set(loci_df.columns)
        if missing:
            logger.error("Loci file missing required columns: %s", ", ".join(sorted(missing)))
            return 2
        if "z_score" in loci_df.columns:
            loci_df = loci_df.sort_values("z_score", ascending=False)
        if args.loci_limit:
            loci_df = loci_df.head(args.loci_limit)
        for _, row in loci_df.iterrows():
            loci.append(
                LocusKey(
                    gene_name=str(row["gene_name"]),
                    match_type=str(row["match_type"]),
                    orf_strand=str(row["orf_strand"]),
                    orf_frame=int(row["orf_frame"]),
                    bin_start=int(row["bin_start"]),
                    bin_end=int(row["bin_end"]),
                )
            )
            frame_values.add(int(row["orf_frame"]))

        if hits_path.exists():
            unique_counts, hit_counts, _ = scan_hits_for_loci(hits_path, args.geom_bins)
    else:
        if not hits_path.exists():
            logger.error("Hits file not found: %s", hits_path)
            return 2

        unique_counts, hit_counts, frame_values = scan_hits_for_loci(
            hits_path,
            args.geom_bins,
        )
        if not unique_counts:
            logger.error("No loci found in hits file.")
            return 2

    if args.frame_mode == "auto":
        if 0 in frame_values:
            frame_mode = "zero"
        elif 3 in frame_values:
            frame_mode = "one"
        else:
            frame_mode = "zero"
            logger.warning("Frame mode auto: no 0/+3 observed; defaulting to zero-based.")
    else:
        frame_mode = args.frame_mode
    logger.info("Using frame mode: %s", frame_mode)

    if not loci:
        loci = select_top_loci(unique_counts, hit_counts, args.min_genomes, args.max_candidates)
    if not loci:
        logger.error("No loci meet min-genomes filter.")
        return 2

    gene_names = {key.gene_name for key in loci}
    logger.info("Selected loci: %d | genes: %d", len(loci), len(gene_names))

    gene_instances = extract_gene_instances(args.gbk_dir, gene_names, args.gene_field)
    if not gene_instances:
        logger.error("No gene instances found for selected loci.")
        return 2

    rng = random.Random(args.seed)

    rows: List[dict] = []
    total_loci = len(loci)
    for idx, locus in enumerate(loci, start=1):
        if idx == 1 or idx % 10 == 0 or idx == total_loci:
            logger.info("Processing locus %d/%d (%s)", idx, total_loci, locus.gene_name)
        instances = gene_instances.get(locus.gene_name, [])
        if not instances:
            continue
        windows = prepare_locus_windows(locus, instances, args.geom_bins)
        if not windows:
            continue
        obs_survival, obs_identity, total, survived, _ = observed_for_locus(
            locus,
            windows,
            frame_mode,
        )
        null_survival, null_identity = null_for_locus(
            locus,
            windows,
            rng,
            args.null_iterations,
            frame_mode,
        )

        null_survival_mean, null_survival_std = summarize_null(null_survival)
        null_identity_mean, null_identity_std = summarize_null(null_identity)

        if null_identity_std and not math.isnan(obs_identity) and null_identity_std > 0:
            z_score = (obs_identity - null_identity_mean) / null_identity_std
        else:
            z_score = float("nan")

        if math.isnan(obs_survival):
            p_survival = float("nan")
        else:
            p_survival = p_value(sum(1 for x in null_survival if x >= obs_survival), len(null_survival))
        if math.isnan(obs_identity):
            p_identity = float("nan")
        else:
            p_identity = p_value(sum(1 for x in null_identity if x >= obs_identity), len(null_identity))

        rows.append({
            "gene_name": locus.gene_name,
            "match_type": locus.match_type,
            "orf_strand": locus.orf_strand,
            "orf_frame": locus.orf_frame,
            "bin_start": locus.bin_start,
            "bin_end": locus.bin_end,
            "genomes": unique_counts.get(locus, 0),
            "hits": hit_counts.get(locus, 0),
            "observed_n": total,
            "observed_survived": survived,
            "obs_survival": obs_survival,
            "obs_identity": obs_identity,
            "null_survival_mean": null_survival_mean,
            "null_survival_std": null_survival_std,
            "null_identity_mean": null_identity_mean,
            "null_identity_std": null_identity_std,
            "z_score": z_score,
            "p_survival": p_survival,
            "p_identity": p_identity,
            "null_samples_survival": len(null_survival),
            "null_samples_identity": len(null_identity),
        })

    if not rows:
        logger.error("No loci scored.")
        return 2

    out_dir = args.output_dir
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / "altframe_constraint_conservation.tsv"
    pd.DataFrame(rows).sort_values(by=["z_score"], ascending=False).to_csv(out_path, sep="\t", index=False)
    logger.info("Wrote %d rows to %s", len(rows), out_path)

    return 0
