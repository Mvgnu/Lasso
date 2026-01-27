from lasso_workbench.core.ranking import TwoStageRanker
from lasso_workbench.schemas.pipeline import ScoredCandidate, RankingConfig


def _cand(candidate_id: str, best: float, top_n_mean: float) -> ScoredCandidate:
    return ScoredCandidate(
        candidate_id=candidate_id,
        record_id="rec1",
        protein_sequence="MKKQTF",
        dna_sequence="ATGAAGAAAACGTTT",
        aa_length=6,
        strand="+",
        genomic_start=0,
        genomic_end=18,
        best_similarity=best,
        best_match_id="ref1",
        top_n_mean_similarity=top_n_mean,
    )


def test_ranking_uses_top_n_mean_without_fallback() -> None:
    config = RankingConfig(score_mode="top_n_mean", min_rule_score=0.0)
    ranker = TwoStageRanker(config)

    cand_a = _cand("a", best=0.9, top_n_mean=0.0)
    cand_b = _cand("b", best=0.1, top_n_mean=0.2)

    ranked = ranker.rank_candidates([cand_a, cand_b])
    assert ranked[0].candidate_id == "b"
    assert ranked[0].embedding_score == 0.2
