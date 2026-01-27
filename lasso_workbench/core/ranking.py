"""
Embedding-first ranking with optional rule-score cutoff.
License: MIT
Author: Magnus Ohle

Stage 1: Sort by embedding score (best_similarity or top_n_mean).
Stage 2: If configured, compute rule scores and drop candidates below a cutoff.

The rule score is used only as a filter (no weighting or normalization).
"""

from typing import List, Optional

import logging

from lasso_workbench.schemas.pipeline import ScoredCandidate, RankingConfig
from lasso_workbench.core.prediction import CorePredictor


logger = logging.getLogger(__name__)


class TwoStageRanker:
    """
    Embedding-first ranking with optional rule-score cutoff.

    - Always sort by embedding score.
    - Optionally drop candidates whose rule score is below a cutoff.
    """

    def __init__(
        self,
        config: Optional[RankingConfig] = None,
        predictor: Optional[CorePredictor] = None,
    ):
        """
        Initialize ranker with configuration.

        Args:
            config: Ranking configuration. Uses defaults if None.
            predictor: Optional CorePredictor for rule-based scoring.
        """
        self.config = config or RankingConfig()
        self._predictor = predictor

    def rank_candidates(
        self,
        candidates: List[ScoredCandidate],
    ) -> List[ScoredCandidate]:
        """
        Rank candidates using embedding order with optional rule cutoff.

        Args:
            candidates: List of scored candidates for a single BGC

        Returns:
            List of ScoredCandidate sorted by final rank (best first)
        """
        if not candidates:
            return []

        # Stage 1: Sort by embedding score (primary signal)
        ordered = sorted(candidates, key=self._embedding_score, reverse=True)

        for cand in ordered:
            score = self._embedding_score(cand)
            cand.embedding_score = score
            cand.combined_score = score
            cand.rule_score_raw = None

        # Stage 2: Optional rule-score cutoff
        if self.config.min_rule_score > 0:
            predictor = self._get_predictor()
            kept: List[ScoredCandidate] = []
            for cand in ordered:
                pred = predictor.predict(
                    cand.protein_sequence,
                    top_n=1,
                    allow_inverted=self.config.allow_inverted_rules,
                )
                raw = float(pred.best_prediction.score) if pred.best_prediction else 0.0
                cand.rule_score_raw = raw
                cand.rule_orientation = pred.best_prediction.orientation if pred.best_prediction else None
                if raw >= self.config.min_rule_score:
                    kept.append(cand)
            ordered = kept

        # Final sort (based on embedding score)
        ordered.sort(key=lambda c: c.embedding_score or 0.0, reverse=True)

        if ordered:
            logger.debug(
                f"Ranked {len(ordered)} candidates: "
                f"top_similarity={self._embedding_score(ordered[0]):.4f}"
            )

        return ordered

    def _embedding_score(self, cand: ScoredCandidate) -> float:
        mode = self.config.score_mode
        if mode == "best_similarity":
            return cand.best_similarity
        if mode == "top_n_mean":
            return cand.top_n_mean_similarity
        # default: return best similarity if no mode is specified
        return cand.best_similarity

    def _get_predictor(self) -> CorePredictor:
        if self._predictor is None:
            self._predictor = CorePredictor()
        return self._predictor
