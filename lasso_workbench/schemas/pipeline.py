from typing import List, Optional, Dict, Any, Literal

from pydantic import BaseModel, Field, field_validator, model_validator, ConfigDict

Strand = Literal["+", "-", "?"]
Device = Literal["auto", "cpu", "cuda", "mps"]
_START_CODONS = {"ATG", "GTG", "TTG"}

class AnnotatedCDS(BaseModel):
    """A CDS feature annotated in a GBK file."""
    # Coordinates are 0-based, half-open: [start, end)
    start: int
    end: int
    strand: Strand
    translation: Optional[str] = None
    protein_id: Optional[str] = None
    locus_tag: Optional[str] = None
    product: Optional[str] = None
    gene_name: Optional[str] = None

    @property
    def length(self) -> int:
        return max(0, self.end - self.start)

class BGCSegment(BaseModel):
    """A BGC segment for analysis."""
    record_id: str
    index: int
    sequence: str  # Storing as str for serialization, convert to Seq when needed
    annotated_cds: List[AnnotatedCDS] = Field(default_factory=list)
    is_lasso: bool = False
    source_file: str = ""
    metadata: Dict[str, Any] = Field(default_factory=dict)

    @property
    def length(self) -> int:
        return len(self.sequence)

class CandidateORF(BaseModel):
    """A candidate ORF extracted from 6-frame translation (Analysis ready)."""
    # Coordinates are 0-based, half-open: [genomic_start, genomic_end)
    candidate_id: str
    record_id: str
    strand: Strand
    frame: int
    aa_length: int
    nt_length: int
    protein_sequence: str
    dna_sequence: str
    genomic_start: int
    genomic_end: int
    start_codon: str = "ATG"
    is_alt_start: bool = False

    @field_validator("start_codon")
    @classmethod
    def _validate_start_codon(cls, value: str) -> str:
        codon = value.upper()
        if codon not in _START_CODONS:
            raise ValueError(f"Invalid start codon: {value}")
        return codon

    @model_validator(mode="after")
    def _validate_lengths(self):
        if self.aa_length != len(self.protein_sequence):
            raise ValueError("aa_length does not match protein_sequence length")
        if self.nt_length != self.aa_length * 3:
            raise ValueError("nt_length does not match aa_length * 3")
        return self

class TransientORF(BaseModel):
    """Raw ORF record from translation (Extraction phase)."""
    # Coordinates are 0-based, half-open: [genomic_start, genomic_end)
    dna: str
    protein: str
    strand: Strand
    frame: int
    genomic_start: int
    genomic_end: int
    aa_len: int
    nt_len: int
    start_codon: str = "ATG"
    
    @property
    def is_alt_start(self) -> bool:
        """Returns True if this ORF uses an alternative start codon (GTG/TTG)."""
        # Note: This duplicates logic but keeps the model self-contained.
        return self.start_codon in _START_CODONS and self.start_codon != "ATG"

    @field_validator("start_codon")
    @classmethod
    def _validate_start_codon(cls, value: str) -> str:
        codon = value.upper()
        if codon not in _START_CODONS:
            raise ValueError(f"Invalid start codon: {value}")
        return codon

class ScoredCandidate(BaseModel):
    """A candidate with ESM-2 similarity score."""
    # Coordinates are 0-based, half-open: [genomic_start, genomic_end)
    candidate_id: str
    record_id: str
    protein_sequence: str
    dna_sequence: str
    aa_length: int
    strand: Strand
    frame: Optional[int] = None
    genomic_start: int
    genomic_end: int
    best_similarity: float
    best_match_id: str
    top_n_mean_similarity: float = 0.0  # Robust ranking metric
    
    # Extra feature columns
    start_codon: str = "ATG"
    is_alt_start: bool = False

    is_known_precursor: bool = False

    # Ranking metadata (populated when RankingConfig used)
    embedding_score: Optional[float] = None  # Score used for stage 1 ranking
    combined_score: float = 0.0  # Final combined score used for ranking
    rule_score_raw: Optional[float] = None  # Raw rule-engine score
    rule_orientation: Optional[Literal["leader_core", "core_leader"]] = None
    genome_count: Optional[int] = None

    @model_validator(mode="after")
    def _validate_lengths(self):
        if self.aa_length != len(self.protein_sequence):
            raise ValueError("aa_length does not match protein_sequence length")
        return self

class RankingConfig(BaseModel):
    """Embedding ordering with optional rule-score cutoff."""
    score_mode: Literal["best_similarity", "top_n_mean"] = Field(
        default="top_n_mean",
        description="Embedding score used for sorting",
    )
    min_rule_score: float = Field(
        default=0.0,
        ge=0.0,
        description="Minimum rule score to keep a candidate (0 disables rule cutoff)",
    )
    allow_inverted_rules: bool = Field(
        default=False,
        description="Allow core_leader (inverted) rule scoring (experimental)",
    )

    model_config = ConfigDict(extra="forbid")


class PipelineMetadata(BaseModel):
    """Provenance and configuration for a pipeline run."""
    model_name: str
    device: Device
    validated_faa: str
    min_aa: int
    max_aa: int
    embed_batch_size: int
    score_batch_size: int
    top_n_mean: int = 5
    ranking_config: Optional["RankingConfig"] = None

    model_config = ConfigDict(extra="forbid")


class PipelineResult(BaseModel):
    """Complete pipeline result for a single BGC."""
    record_id: str
    is_lasso: bool
    segment_length_nt: int
    num_annotated_cds: int
    num_candidates: int
    candidates: List[ScoredCandidate]
    annotated_cds: List[AnnotatedCDS]
    source_file: str
    pipeline_metadata: PipelineMetadata
