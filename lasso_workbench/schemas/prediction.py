from typing import List, Dict, Optional, Literal

from pydantic import BaseModel, Field, ConfigDict

ConfidenceLevel = Literal["low", "medium", "high"]
Orientation = Literal["leader_core", "core_leader"]

class PredictionResult(BaseModel):
    """Result of core prediction for a single cleavage site."""
    cleavage_site: int
    leader: str
    core: str
    score: float
    score_breakdown: Dict[str, float] = Field(default_factory=dict)
    confidence: ConfidenceLevel = "medium"
    reasons: List[str] = Field(default_factory=list) #TODO: currently not appending rule/reason, missing connector
    orientation: Orientation = "leader_core"

    model_config = ConfigDict(frozen=True)

class PrecursorPrediction(BaseModel):
    """Full prediction result for a precursor sequence."""
    sequence: str
    length: int
    predictions: List[PredictionResult]
    best_prediction: Optional[PredictionResult] = None
    
    def __post_init__(self):
        # Pydantic doesn't use __post_init__ the same way dataclasses do, 
        # but we can use a root validator or model_validator if we want computed fields.
        # For now, we'll leave it as a data structure and let the service populate it.
        pass

class PredictionRequest(BaseModel):
    """Request object for prediction."""
    sequence: str
    top_n: int = 3
    min_leader: int = 10
    min_core: Optional[int] = None

    
