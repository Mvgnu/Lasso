"""
Lasso peptide discovery pipeline modules.

This package contains the full pipeline for:
1. ORF harvesting from GenBank/antiSMASH files via 6-frame translation
2. ESM-2 embedding of precursor candidates  
3. Semantic scoring and ranking against validated precursors
4. Core peptide prediction using rule-based approaches
5. Mass spec target generation for lab validation

Pipeline workflow:
    GBK/FASTA → 6-frame ORF extraction → ESM-2 embedding → 
    → Cosine similarity scoring → Ranking → Core prediction → MS targets

Modules:
- semantic_pipeline: Main unified pipeline for precursor discovery
- orf_extraction: Extract short ORFs via 6-frame translation
- embed_precursor_candidates: ESM-2 embedding and similarity scoring
- generate_mass_spec_targets: Mass spec validation target generation
"""

from lasso_workbench.core.prediction import (
    CorePredictor,
    RuleEngine,
    PrecursorPrediction as CorePrediction,  # Map to similar name
)

# batch module archived
from lasso_workbench.schemas.prediction import PredictionResult as CleavageSite

from lasso_workbench.schemas.pipeline import TransientORF
from lasso_workbench.pipeline.orf_extraction import (
    chunk_orfs,
)

from lasso_workbench.pipeline.semantic_pipeline import (
    run_semantic_pipeline,
    parse_gbk_file,
    parse_gbk_folder,
    extract_orfs,
    BGCSegment,
    CandidateORF,
    ScoredCandidate,
    PipelineResult,
)

__all__ = [
    # Core prediction
    "CorePredictor",
    "RuleEngine",
    "CleavageSite",
    "CorePrediction",
    # BatchCorePredictor and PeptidaseGrouper archived
    # ORF harvesting
    "TransientORF",
    "chunk_orfs",
    # Semantic pipeline
    "run_semantic_pipeline",
    "parse_gbk_file",
    "parse_gbk_folder",
    "extract_orfs",
    "BGCSegment",
    "CandidateORF",
    "ScoredCandidate",
    "PipelineResult",
]
