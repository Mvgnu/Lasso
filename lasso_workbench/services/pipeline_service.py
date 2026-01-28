from pathlib import Path
from typing import List, Tuple, Optional, Dict
import pandas as pd
import tempfile

from lasso_workbench.pipeline import semantic_pipeline as pipeline
from lasso_workbench.pipeline.semantic_pipeline import (
    BGCSegment,
    PipelineResult,
    parse_gbk_file,
    Device,
)
from lasso_workbench.schemas.pipeline import RankingConfig
from lasso_workbench.core.prediction import CorePredictor

class PipelineService:
    """
    Service for running the semantic precursor discovery pipeline.
    """
    def __init__(self):
        pass

    def parse_uploaded_gbk_files(self, file_paths: List[str]) -> Tuple[List[BGCSegment], pd.DataFrame]:
        """Parse multiple GBK files and return segments and summary."""
        segments = []
        summary_rows = []
        
        for idx, path_str in enumerate(file_paths, 1):
            path = Path(path_str)
            # Parse once and pass BGCSegment through preview + run to avoid re-reading GBKs.
            seg = parse_gbk_file(path, idx)
            if seg:
                segments.append(seg)
                summary_rows.append({
                    "File": path.name,
                    "Record ID": seg.record_id,
                    "Length (nt)": seg.length,
                    "CDS Count": len(seg.annotated_cds),
                    "Lasso Detected": "Yes" if seg.is_lasso else "No"
                })
        
        return segments, pd.DataFrame(summary_rows)

    def run_pipeline(
        self,
        segments: List[BGCSegment],
        validated_faa_path: Path,
        ranking_config: RankingConfig,
        min_aa: int = 20,
        max_aa: int = 120,
        device: Device = "auto",
        model_name: Optional[str] = None,
        out_dir: Optional[Path] = None,
        progress_callback=None,
        ranker_predictor: Optional[CorePredictor] = None,
        filter_lasso_only: bool = False,
    ) -> Tuple[List[PipelineResult], pd.DataFrame]:
        """
        Run the full pipeline: ORF extraction -> Embedding -> Scoring -> Reporting.

        Args:
            segments: List of parsed BGC segments
            validated_faa_path: Path to validated precursor FASTA
            min_aa: Minimum ORF length in amino acids
            max_aa: Maximum ORF length in amino acids
            device: Compute device (auto, cpu, cuda, mps)
            model_name: ESM-2 model name
            out_dir: Output directory for intermediate files
            progress_callback: Optional progress callback function
            ranking_config: Two-stage ranking configuration
            ranker_predictor: Optional CorePredictor for rule-based re-ranking
            filter_lasso_only: If True, only analyze lasso-annotated BGCs
        """
        if out_dir is None:
            out_dir = Path(tempfile.mkdtemp(prefix="lasso_pipeline_"))

        results, _ = pipeline.run_semantic_pipeline(
            gbk_files=None,
            segments=segments,
            validated_faa=validated_faa_path,
            output_dir=out_dir,
            min_aa=min_aa,
            max_aa=max_aa,
            device=device,
            model_name=model_name,
            progress=progress_callback,
            embed_batch_size=100,
            score_batch_size=100,
            top_n_mean=5,
            ranking_config=ranking_config,
            ranker_predictor=ranker_predictor,
            filter_lasso_only=filter_lasso_only,
        )
        
        # Create flat dataframe for display
        rows = []
        for res in results:
            for cand in res.candidates:
                score = cand.embedding_score if cand.embedding_score is not None else cand.best_similarity
                rows.append({
                    "Record": cand.record_id,
                    "Candidate": cand.candidate_id,
                    "Score": score,
                    "Rule Score": cand.rule_score_raw,
                    "Sequence": cand.protein_sequence,
                    "DNA Sequence (excl. stop)": cand.dna_sequence,
                    "Strand": cand.strand,
                    # Display 1-based inclusive coordinates for UI readability and consistency with GBK.
                    "Start": cand.genomic_start + 1,
                    "End": cand.genomic_end +1,
                })
        
        return results, pd.DataFrame(rows)
