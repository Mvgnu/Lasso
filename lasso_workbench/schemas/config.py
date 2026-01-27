from typing import Dict, List, Optional, Any, Type
from pydantic import BaseModel, Field, ConfigDict, field_validator

class CoreStartResidueParameters(BaseModel):
    scores: Dict[str, float]
    penalty_other: float = -0.5
    hard_forbidden: List[str] = Field(default_factory=list)
    penalty_forbidden: float = -2.0

    model_config = ConfigDict(extra="forbid")


class RingClosureResidueParameters(BaseModel):
    allowed_positions: List[int]
    residues: List[str]
    position_weights: Dict[str, float] = Field(default_factory=dict)
    penalty_if_missing_in_allowed_window: float = -5.0

    model_config = ConfigDict(extra="forbid")


class CoreLengthParameters(BaseModel):
    min: int
    max: int
    optimal_min: int
    optimal_max: int
    score_optimal: float = 1.0
    score_acceptable: float = 0.5
    penalty_outside: float = -1.0

    model_config = ConfigDict(extra="forbid")


class LeaderLengthParameters(BaseModel):
    min: int
    max: int
    optimal_min: Optional[int] = None
    optimal_max: Optional[int] = None

    model_config = ConfigDict(extra="forbid")


class CleavageMotifParameters(BaseModel):
    motifs: Dict[str, float] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")


class LeaderMotifsParameters(BaseModel):
    motifs: List[str]

    model_config = ConfigDict(extra="forbid")


class StericLockRegionParameters(BaseModel):
    bulky_residues: List[str]
    relative_to_acceptor: Dict[str, int]
    require_at_least: int = 1

    model_config = ConfigDict(extra="forbid")


class LeaderPenultimateResidueParameters(BaseModel):
    preferred: Dict[str, float] = Field(default_factory=dict)
    allowed_rare: Dict[str, float] = Field(default_factory=dict)
    penalty_other: float = -0.5

    model_config = ConfigDict(extra="forbid")


class LeaderSignatureMotifsParameters(BaseModel):
    patterns: Dict[str, float] = Field(default_factory=dict)
    allow_regex: bool = False

    model_config = ConfigDict(extra="forbid")


class PtmContextHintsParameters(BaseModel):
    if_kinase_present: Dict[str, List[str]] = Field(default_factory=dict)
    if_glycosyltransferase_present: Dict[str, List[str]] = Field(default_factory=dict)
    if_acyltransferase_present: Dict[str, List[str]] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")


class GlycineRichCoreParameters(BaseModel):
    min_ratio: float

    model_config = ConfigDict(extra="forbid")


class HydrophobicityTransitionParameters(BaseModel):
    min_difference: float
    scale: Optional[str] = None
    window_size: Optional[int] = None

    model_config = ConfigDict(extra="forbid")


class ChargeTransitionParameters(BaseModel):
    leader_more_charged: bool
    min_difference: float = 0.1

    model_config = ConfigDict(extra="forbid")


class TailLockResidueParameters(BaseModel):
    region_start: int
    residues: List[str]
    require_at_least: int = 1

    model_config = ConfigDict(extra="forbid")


class ForbiddenCorePatternsParameters(BaseModel):
    patterns: List[str]
    penalty_per_match: float = -1.0

    model_config = ConfigDict(extra="forbid")


class AlignmentConservationParameters(BaseModel):
    min_sequences: int
    conservation_threshold: float
    method: str

    model_config = ConfigDict(extra="forbid")


class Position2ConstraintParameters(BaseModel):
    preferred: List[str]
    forbidden: List[str]
    neutral: List[str]

    model_config = ConfigDict(extra="forbid")


RuleParameters = (
    CoreStartResidueParameters
    | RingClosureResidueParameters
    | CoreLengthParameters
    | LeaderLengthParameters
    | CleavageMotifParameters
    | LeaderMotifsParameters
    | StericLockRegionParameters
    | LeaderPenultimateResidueParameters
    | LeaderSignatureMotifsParameters
    | GlycineRichCoreParameters
    | HydrophobicityTransitionParameters
    | ChargeTransitionParameters
    | TailLockResidueParameters
    | ForbiddenCorePatternsParameters
    | AlignmentConservationParameters
    | Position2ConstraintParameters
    | PtmContextHintsParameters
)


_PARAMETER_MODELS: Dict[str, Type[BaseModel]] = {
    "core_start_residue": CoreStartResidueParameters,
    "ring_closure_residue": RingClosureResidueParameters,
    "core_length": CoreLengthParameters,
    "leader_length": LeaderLengthParameters,
    "cleavage_motif": CleavageMotifParameters,
    "leader_motifs": LeaderMotifsParameters,
    "steric_lock_region": StericLockRegionParameters,
    "leader_penultimate_residue": LeaderPenultimateResidueParameters,
    "leader_signature_motifs": LeaderSignatureMotifsParameters,
    "glycine_rich_core": GlycineRichCoreParameters,
    "hydrophobicity_transition": HydrophobicityTransitionParameters,
    "charge_transition": ChargeTransitionParameters,
    "tail_lock_residue": TailLockResidueParameters,
    "forbidden_core_patterns": ForbiddenCorePatternsParameters,
    "alignment_conservation": AlignmentConservationParameters,
    "position_2_constraint": Position2ConstraintParameters,
    "ptm_context_hints": PtmContextHintsParameters,
}


class RuleParameter(BaseModel):
    """Configuration for a specific rule parameter."""
    name: str
    description: Optional[str] = None
    sources: List[str] = Field(default_factory=list)
    enabled: bool = True
    weight: float = 1.0
    parameters: RuleParameters
    ui: Dict[str, Any] = Field(default_factory=dict)
    # UI hints can be stored here or separate, keeping it simple for now

    model_config = ConfigDict(extra="forbid")

    @field_validator("parameters", mode="before")
    @classmethod
    def _validate_parameters(cls, value: Any, info) -> RuleParameters:
        rule_name = info.data.get("name")
        if not rule_name:
            raise ValueError("Rule name is required for parameter validation")
        model = _PARAMETER_MODELS.get(rule_name)
        if model is None:
            raise ValueError(f"Unknown rule name for parameters: {rule_name}")
        return model.model_validate(value)

class RuleSet(BaseModel):
    """Collection of rules for prediction."""
    version: str = "1.0"
    description: str = ""
    rules: Dict[str, RuleParameter] = Field(default_factory=dict)

    model_config = ConfigDict(extra="forbid")
