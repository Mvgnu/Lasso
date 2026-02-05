#!/usr/bin/env python3
"""
Unified Semantic Precursor Discovery Pipeline

This module provides the core pipeline for lasso peptide precursor discovery:
1. Load nucleotide sequences from GBK files (MiBIG, antiSMASH, or custom)
2. Perform 6-frame translation to extract ALL ORFs (via orf_extraction)
3. Embed candidates with ESM-2 (via embed_precursor_candidates)
4. Rank against validated precursor reference set
5. Generate reports with ALL scores (not just top N)

Designed for GUI integration - works with uploaded files and configurable parameters.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Callable, Set, Literal

import pandas as pd
from Bio.Seq import Seq

from lasso_workbench.schemas.pipeline import (
    AnnotatedCDS,
    BGCSegment,
    CandidateORF,
    ScoredCandidate,
    PipelineResult,
    PipelineMetadata,
    RankingConfig,
)
from lasso_workbench.core.ranking import TwoStageRanker
from lasso_workbench.core.prediction import CorePredictor
from lasso_workbench.pipeline.orf_extraction import chunk_orfs
from lasso_workbench.pipeline.embed_precursor_candidates import (
    score_candidates_in_memory,
    load_fasta,
)
from lasso_workbench.pipeline.esm_embedder import ESM2Embedder
from lasso_workbench.utils.sequence_io import read_genbank_record, write_fasta_pairs

logger = logging.getLogger(__name__)

Device = Literal["auto", "cpu", "cuda", "mps"]


@dataclass(frozen=True)
class PipelineRunConfig:
    min_aa: int = 20
    max_aa: int = 120
    device: Device = "auto"
    model_name: Optional[str] = None
    embed_batch_size: int = 100
    score_batch_size: int = 100
    top_n_mean: int = 5


@dataclass(frozen=True)
class PipelineContext:
    validated_faa: Path
    output_dir: Path
    segments: List["BGCSegment"]


# ---------------------------------------------------------------------------
# GBK Parsing
# ---------------------------------------------------------------------------

def parse_gbk_file(gbk_path: Path, index: int = 1) -> Optional[BGCSegment]:
    """
    Parse a GenBank file and extract nucleotide sequence + CDS annotations.
    
    Works with MiBIG, antiSMASH, or any standard GBK file.
    """
    record = read_genbank_record(gbk_path)
    if record is None:
        logger.warning(f"Could not parse {gbk_path}")
        return None

    sequence = record.seq

    record_id = record.id or gbk_path.stem
    
    # Extract annotated CDS
    annotated_cds: List[AnnotatedCDS] = []
    is_lasso = False
    region_ranges: List[Tuple[int, int]] = []

    for feat in record.features:
        if feat.type == "CDS":
            translation = feat.qualifiers.get("translation", [None])[0]
            protein_id = feat.qualifiers.get("protein_id", [None])[0]
            locus_tag = feat.qualifiers.get("locus_tag", [None])[0]
            product = feat.qualifiers.get("product", [None])[0]
            gene_name = feat.qualifiers.get("gene", [None])[0]
            if feat.location.strand == 1:
                strand = "+"
            elif feat.location.strand == -1:
                strand = "-"
            else:
                strand = "?"

            annotated_cds.append(AnnotatedCDS(
                start=int(feat.location.start),
                end=int(feat.location.end),
                strand=strand,
                translation=translation,
                protein_id=protein_id,
                locus_tag=locus_tag,
                product=product,
                gene_name=gene_name,
            ))
        
        if feat.type in {"region", "cluster"}:
            products = feat.qualifiers.get("product", [])
            if not products:
                products = feat.qualifiers.get("products", [])
            if any("lasso" in str(p).lower() for p in products):
                start = int(feat.location.start)
                end = int(feat.location.end)
                if start < end:
                    region_ranges.append((start, end))
                    is_lasso = True

    # BGCSegment is the shared data object between UI preview and pipeline run.
    # It carries parsed sequence/annotations so we don't re-parse GBKs twice.
    return BGCSegment(
        record_id=record_id,
        index=index,
        sequence=str(sequence),
        annotated_cds=annotated_cds,
        is_lasso=is_lasso,
        source_file=str(gbk_path),
        metadata={
            "length_nt": len(sequence),
            "num_annotated_cds": len(annotated_cds),
            **(
                {
                    "region_ranges": [{"start": s, "end": e} for s, e in region_ranges],
                    "region_source": "lasso",
                }
                if region_ranges
                else {}
            ),
        },
    )


def parse_gbk_folder(folder_path: Path) -> List[BGCSegment]:
    """Parse all GBK files in a folder."""
    segments: List[BGCSegment] = []
    
    for idx, gbk_file in enumerate(sorted(folder_path.glob("*.gbk")), start=1):
        segment = parse_gbk_file(gbk_file, idx)
        if segment:
            segments.append(segment)

    return ensure_unique_record_ids(segments)


# ---------------------------------------------------------------------------
# ORF Extraction (6-frame translation)
# ---------------------------------------------------------------------------

def extract_orfs(
    segments: List[BGCSegment],
    min_aa: int = 20,
    max_aa: int = 120,
) -> Tuple[List[CandidateORF], List[dict]]:
    """
    Extract ORFs from all segments and return candidates + metadata rows.
    """
    all_candidates: List[CandidateORF] = []
    metadata_rows: List[dict] = []

    for segment in segments:
        seq = Seq(segment.sequence)
        region_ranges = segment.metadata.get("region_ranges", [])
        use_region_ranges = bool(region_ranges)
        if segment.metadata.get("region_source") == "lasso" and not segment.metadata.get("restrict_to_lasso", False):
            use_region_ranges = False
        if use_region_ranges:
            windows = [(int(item["start"]), int(item["end"])) for item in region_ranges]
        else:
            windows = [(0, len(seq))]
        candidate_index = 1

        for window_start, window_end in windows:
            window_seq = seq[window_start:window_end]
            # Forward strand orfs (frames 0/1/2)
            forward_orfs = list(chunk_orfs(window_seq, "+", window_start, window_end, min_aa, max_aa))

            # Reverse-complement strand orfs (frames 0/1/2), pass "-" to correctly map back to GBK
            rc_seq = window_seq.reverse_complement()
            reverse_orfs = list(chunk_orfs(rc_seq, "-", window_start, window_end, min_aa, max_aa))

            orfs = forward_orfs + reverse_orfs
            for orf in orfs:
                candidate_id = f"{segment.record_id}|orf{candidate_index:04d}"
                candidate_index += 1
                cand = CandidateORF(
                    candidate_id=candidate_id,
                    record_id=segment.record_id,
                    strand=orf.strand,
                    frame=int(orf.frame),
                    aa_length=int(orf.aa_len),
                    nt_length=int(orf.nt_len),
                    protein_sequence=orf.protein,
                    dna_sequence=orf.dna, 
                    genomic_start=int(orf.genomic_start),
                    genomic_end=int(orf.genomic_end),
                    start_codon=str(getattr(orf, "start_codon", "ATG")),
                    is_alt_start=bool(getattr(orf, "is_alt_start", False)),
                )
                all_candidates.append(cand)
                metadata_rows.append({
                    "candidate_id": cand.candidate_id,
                    "record_id": cand.record_id,
                    "strand": cand.strand,
                    "frame": cand.frame,
                    "aa_length": cand.aa_length,
                    "protein_sequence": cand.protein_sequence,
                    "dna_sequence": cand.dna_sequence,
                    "genomic_start": cand.genomic_start,
                    "genomic_end": cand.genomic_end,
                })

    return all_candidates, metadata_rows

def write_fasta(records: List[Tuple[str, str]], path: Path) -> None:
    """Write sequences to FASTA format using BioPython."""
    write_fasta_pairs(records, path)


def candidates_to_fasta(candidates: List[CandidateORF]) -> List[Tuple[str, str]]:
    """Convert candidates to FASTA (id, sequence) pairs."""
    return [(c.candidate_id, c.protein_sequence) for c in candidates]


def build_pipeline_metadata(
    model_name: str,
    device: Device,
    validated_faa: Path,
    min_aa: int,
    max_aa: int,
    embed_batch_size: int,
    score_batch_size: int,
    top_n_mean: int,
    ranking_config: RankingConfig,
) -> PipelineMetadata:
    """Build PipelineMetadata with standardized fields."""
    return PipelineMetadata(
        model_name=model_name or "",
        device=device,
        validated_faa=str(validated_faa),
        min_aa=min_aa,
        max_aa=max_aa,
        embed_batch_size=embed_batch_size,
        score_batch_size=score_batch_size,
        top_n_mean=top_n_mean,
        ranking_config=ranking_config,
    )


def ensure_unique_record_ids(
    segments: List[BGCSegment],
    progress: Optional[Callable[[str], None]] = None,
) -> List[BGCSegment]:
    """
    Normalize and de-duplicate record IDs across segments.

    Duplicate record IDs can cause candidate-id collisions and result-map overwrites.
    """
    seen: Dict[str, int] = {}
    for idx, segment in enumerate(segments, start=1):
        base_id = (segment.record_id or "").strip() or f"segment_{idx}"
        count = seen.get(base_id, 0) + 1
        seen[base_id] = count

        if count == 1:
            segment.record_id = base_id
            continue

        new_id = f"{base_id}__{count}"
        segment.record_id = new_id
        if progress is not None:
            progress(f"Duplicate record_id '{base_id}' renamed to '{new_id}'")
    return segments


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Report Building
# ---------------------------------------------------------------------------

def build_pipeline_results(
    similarity_rows: List[dict],
    metadata_rows: List[dict],
    segments: Dict[str, BGCSegment],
    ranking_config: RankingConfig,
    known_precursors: Optional[List[str]] = None,
    save_all_scores: bool = True,
    pipeline_metadata: Optional[PipelineMetadata] = None,
    ranker_predictor: Optional[CorePredictor] = None,
) -> List[PipelineResult]:
    """
    Build structured results for each BGC.
    """
    if known_precursors:
        known_upper = [seq.upper() for seq in known_precursors if seq]
    else:
        known_upper = []

    def _is_known_precursor(sequence: str) -> bool:
        if not known_upper:
            return False
        seq_upper = sequence.upper()
        return any(precursor in seq_upper for precursor in known_upper)

    metadata_by_id = {row["candidate_id"]: row for row in metadata_rows}
    merged_by_record: Dict[str, List[dict]] = {}
    for row in similarity_rows:
        candidate_id = row.get("candidate_id")
        # some sanity checks
        if not candidate_id:
            raise ValueError("Missing candidate_id in similarity rows")
        meta = metadata_by_id.get(candidate_id)
        if not meta:
            raise ValueError(f"Missing metadata for candidate_id={candidate_id}")
        merged = {**meta, **row}
        record_id = merged.get("record_id")
        if not record_id:
            raise ValueError(f"Missing record_id for candidate_id={candidate_id}")
        merged_by_record.setdefault(record_id, []).append(merged)

    results: List[PipelineResult] = []

    for record_id, group in merged_by_record.items():
        segment = segments.get(record_id)
        if not segment:
            continue
        
        # Sort by similarity score descending (default to top_n_mean)
        if ranking_config.score_mode == "best_similarity":
            score_key = "best_similarity"
        elif group and "top_n_mean_similarity" in group[0]:
            score_key = "top_n_mean_similarity"
        else:
            score_key = "best_similarity"

        if any(score_key not in row for row in group):
            raise ValueError(f"Missing {score_key} in similarity rows for record_id={record_id}")
        sorted_group = sorted(group, key=lambda row: float(row[score_key]), reverse=True)
        
        candidates: List[ScoredCandidate] = []
        for row in sorted_group:
            candidates.append(ScoredCandidate(
                candidate_id=row["candidate_id"],
                record_id=record_id,
                protein_sequence=row["protein_sequence"],
                dna_sequence=row["dna_sequence"],
                aa_length=int(row["aa_length"]),
                strand=row["strand"],
                frame=int(row["frame"]) if row.get("frame") is not None else None,
                genomic_start=int(row["genomic_start"]),
                genomic_end=int(row["genomic_end"]),
                best_similarity=float(row["best_similarity"]),
                best_match_id=row["best_match_id"],
                is_known_precursor=_is_known_precursor(row["protein_sequence"]),
                top_n_mean_similarity=float(row.get("top_n_mean_similarity", 0.0)),
            ))

        # Apply two-stage ranking (embedding order + optional rule cutoff)
        if candidates:
            ranker = TwoStageRanker(ranking_config, predictor=ranker_predictor)
            candidates = ranker.rank_candidates(candidates)
            logger.debug(
                f"Applied two-stage ranking to {record_id}: "
                f"{len(candidates)} candidates"
            )

        if pipeline_metadata is None:
            raise ValueError("pipeline_metadata is required")
        results.append(PipelineResult(
            record_id=record_id,
            is_lasso=segment.is_lasso,
            segment_length_nt=segment.length,
            num_annotated_cds=len(segment.annotated_cds),
            num_candidates=len(candidates),
            candidates=candidates,
            annotated_cds=segment.annotated_cds,
            source_file=segment.source_file,
            pipeline_metadata=pipeline_metadata,
        ))
    
    return results


def results_to_json(results: List[PipelineResult]) -> List[Dict]:
    """Convert results to JSON-serializable format."""
    output = []

    cand_fields = {
        "candidate_id",
        "protein_sequence",
        "dna_sequence",
        "aa_length",
        "strand",
        "frame",
        "genomic_start",
        "genomic_end",
        "best_similarity",
        "best_match_id",
        "top_n_mean_similarity",
        "embedding_score",
        "is_known_precursor",
        "rule_orientation",
        "rule_score_raw",
        "genome_count",
    }

    cds_fields = {
        "start",
        "end",
        "strand",
        "translation",
        "protein_id",
        "locus_tag",
        "product",
    }

    for result in results:
        candidates_json = []
        for cand in result.candidates:
            cand_data = cand.model_dump()
            candidate = {key: cand_data.get(key) for key in cand_fields}
            candidates_json.append(candidate)

        cds_json = []
        for cds in result.annotated_cds:
            cds_data = cds.model_dump()
            cds_entry = {key: cds_data.get(key) for key in cds_fields}
            cds_entry["length_nt"] = cds.length
            cds_json.append(cds_entry)

        output.append({
            "record_id": result.record_id,
            "is_lasso": result.is_lasso,
            "segment_length_nt": result.segment_length_nt,
            "num_annotated_cds": result.num_annotated_cds,
            "num_candidates": result.num_candidates,
            "candidates": candidates_json,
            "annotated_cds": cds_json,
            "source_file": result.source_file,
            "coordinate_system": "0-based half-open",
            "pipeline_metadata": result.pipeline_metadata.model_dump(),
        })

    return output


def results_to_dataframe(results: List[PipelineResult]) -> pd.DataFrame:
    """Flatten results to a single DataFrame with candidates."""
    frames: List[pd.DataFrame] = []

    for result in results:
        if not result.candidates:
            continue

        cand_df = pd.DataFrame([cand.model_dump() for cand in result.candidates])
        cand_df["is_lasso_bgc"] = result.is_lasso

        frames.append(cand_df)

    if not frames:
        return pd.DataFrame()

    df = pd.concat(frames, ignore_index=True)
    preferred_columns = [
        "record_id",
        "is_lasso_bgc",
        "candidate_id",
        "protein_sequence",
        "dna_sequence",
        "aa_length",
        "frame",
        "strand",
        "genomic_start",
        "genomic_end",
        "start_codon",
        "is_alt_start",
        "best_similarity",
        "best_match_id",
        "top_n_mean_similarity",
        "embedding_score",
        "combined_score",
        "rule_score_raw",
        "rule_orientation",
        "is_known_precursor",
        "genome_count",
    ]
    ordered = [c for c in preferred_columns if c in df.columns]
    remaining = [c for c in df.columns if c not in ordered]
    return df.loc[:, ordered + remaining]


# ---------------------------------------------------------------------------
# Main Pipeline Function
# ---------------------------------------------------------------------------

def run_semantic_pipeline(
    gbk_files: Optional[List[Path]],
    validated_faa: Path,
    output_dir: Path,
    ranking_config: RankingConfig,
    min_aa: int = 20,
    max_aa: int = 120,
    device: Optional[str] = None,
    model_name: Optional[str] = None,
    progress: Optional[Callable[[str], None]] = None,
    embed_batch_size: int = 100,
    score_batch_size: int = 100,
    top_n_mean: int = 5,
    ranker_predictor: Optional[CorePredictor] = None,
    segments: Optional[List[BGCSegment]] = None,
    filter_lasso_only: bool = False,
    orf_index_db: Optional[Path] = None,
) -> Tuple[List[PipelineResult], pd.DataFrame]:
    """
    Run the full semantic precursor discovery pipeline.
    """

    if not model_name:
        raise ValueError("model_name must be provided")
    
    config = PipelineRunConfig(
        min_aa=min_aa,
        max_aa=max_aa,
        device=(device or "auto"),
        model_name=model_name,
        embed_batch_size=embed_batch_size,
        score_batch_size=score_batch_size,
        top_n_mean=top_n_mean,
    )

    emit = progress if progress is not None else logger.info

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Step 1: Parse GBK files (or use provided segments), inform UI
    emit(f"Embedding model: {config.model_name}")
    emit(f"Embedding device requested: {config.device}")
    emit(f"Embedding batch_size: {config.embed_batch_size} | scoring chunk: {config.score_batch_size}")
    
    segments_list: List[BGCSegment] = []
    if segments is not None:
        segments_list = list(segments)
        emit(f"[1/6] Using {len(segments_list)} pre-parsed BGC segments")
    else:
        # Handle usage as script (GBK provided as args)
        gbk_files = gbk_files or []
        emit(f"[1/6] Parsing {len(gbk_files)} GBK files")
        for idx, gbk_path in enumerate(gbk_files, start=1):
            segment = parse_gbk_file(Path(gbk_path), idx)
            if segment:
                segments_list.append(segment)

    if filter_lasso_only:
        before_count = len(segments_list)
        segments_list = [seg for seg in segments_list if seg.is_lasso]
        emit(f"[1/6] Filtered to lasso BGCs: {len(segments_list)} of {before_count}")
        if not segments_list:
            raise ValueError("No lasso-annotated BGCs found after filtering")

    ensure_unique_record_ids(segments_list, progress=emit)

    for seg in segments_list:
        seg.metadata["restrict_to_lasso"] = bool(filter_lasso_only)

    context = PipelineContext(validated_faa=validated_faa, output_dir=output_dir, segments=segments_list)
    
    # Step 2: Extract ORFs OR via 6-frame translation, memory bottleneck atm when a lot of GBKs are provided
    emit(f"[2/6] Harvesting ORFs via 6-frame translation ({config.min_aa}-{config.max_aa} aa)")
    candidates, metadata_rows = extract_orfs(
        segments_list,
        config.min_aa,
        config.max_aa,
    )
    
    if not candidates:
        raise ValueError("No ORFs extracted from BGC sequences")
    
    emit(f"Extracted {len(candidates)} candidate ORFs")
    # Step 3: Embed + score in memory
    emit("[3/6] Embedding candidates + scoring vs validated set")
    validated = load_fasta(context.validated_faa)
    known_precursors = [seq for _, seq in validated]
    candidate_records = candidates_to_fasta(candidates)
    similarity_rows = score_candidates_in_memory(
        validated_sequences=validated,
        candidate_sequences=candidate_records,
        embedder=ESM2Embedder(config.model_name, device=config.device if config.device != "auto" else None),
        score_batch_size=config.score_batch_size,
        embed_batch_size=config.embed_batch_size,
        top_n_mean=config.top_n_mean,
        progress_callback=progress,
    )

    emit("Embedding/scoring done")
    
    # Step 4: Build results
    emit("[4/6] Building per-BGC result objects")

    pipeline_metadata = build_pipeline_metadata(
        model_name=config.model_name or "",
        device=config.device,
        validated_faa=context.validated_faa,
        min_aa=config.min_aa,
        max_aa=config.max_aa,
        embed_batch_size=config.embed_batch_size,
        score_batch_size=config.score_batch_size,
        top_n_mean=config.top_n_mean,
        ranking_config=ranking_config,
    )

    segment_map = {s.record_id: s for s in segments_list}
    results = build_pipeline_results(
        similarity_rows,
        metadata_rows,
        segment_map,
        pipeline_metadata=pipeline_metadata,
        ranking_config=ranking_config,
        ranker_predictor=ranker_predictor,
        known_precursors=known_precursors,
    )

    emit("Built results")

    db_path = orf_index_db or Path("results/lasso_orf_index.sqlite")
    if db_path.exists():
        from lasso_workbench.altframe.indexer import lookup_genome_counts

        sequences = {cand.protein_sequence for res in results for cand in res.candidates}
        genome_counts = lookup_genome_counts(db_path, list(sequences))
        for res in results:
            for cand in res.candidates:
                cand.genome_count = genome_counts.get(cand.protein_sequence, 0)
        emit(f"Attached genome counts from {db_path}")
    else:
        emit(f"No ORF index found at {db_path}; skipping genome count annotation")
    
    # Step 5: Export results
    emit("[5/6] Exporting TSV/JSON outputs")
    all_candidates_df = results_to_dataframe(results)

    emit("Export done")
    emit("[6/6] Pipeline complete")
    emit(f"Pipeline complete. {len(all_candidates_df)} total scored candidates.")
    
    return results, all_candidates_df
