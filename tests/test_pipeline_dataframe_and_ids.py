from lasso_workbench.pipeline.semantic_pipeline import (
    ensure_unique_record_ids,
    results_to_dataframe,
)
from lasso_workbench.schemas.pipeline import (
    BGCSegment,
    PipelineMetadata,
    PipelineResult,
    RankingConfig,
    ScoredCandidate,
)


def test_ensure_unique_record_ids_renames_duplicates():
    seg_a = BGCSegment(record_id="rec", index=1, sequence="ATGC")
    seg_b = BGCSegment(record_id="rec", index=2, sequence="ATGC")
    seg_c = BGCSegment(record_id="rec", index=3, sequence="ATGC")

    ensure_unique_record_ids([seg_a, seg_b, seg_c])

    assert seg_a.record_id == "rec"
    assert seg_b.record_id == "rec__2"
    assert seg_c.record_id == "rec__3"


def test_results_to_dataframe_preserves_extended_candidate_fields():
    candidate = ScoredCandidate(
        candidate_id="rec|orf0001",
        record_id="rec",
        protein_sequence="MKT",
        dna_sequence="ATGAAAACC",
        aa_length=3,
        strand="+",
        frame=0,
        genomic_start=0,
        genomic_end=9,
        best_similarity=0.9,
        best_match_id="ref|locus=x",
        top_n_mean_similarity=0.85,
        embedding_score=0.85,
        combined_score=0.85,
        rule_score_raw=4.2,
        rule_orientation="leader_core",
    )
    result = PipelineResult(
        record_id="rec",
        is_lasso=True,
        segment_length_nt=100,
        num_annotated_cds=0,
        num_candidates=1,
        candidates=[candidate],
        annotated_cds=[],
        source_file="rec.gbk",
        pipeline_metadata=PipelineMetadata(
            model_name="facebook/esm2_t6_8M_UR50D",
            device="cpu",
            validated_faa="validated.faa",
            min_aa=20,
            max_aa=120,
            embed_batch_size=100,
            score_batch_size=100,
            top_n_mean=5,
            ranking_config=RankingConfig(),
        ),
    )

    df = results_to_dataframe([result])

    assert "candidate_id" in df.columns
    assert "combined_score" in df.columns
    assert "rule_score_raw" in df.columns
    assert df.loc[0, "candidate_id"] == "rec|orf0001"
    assert float(df.loc[0, "rule_score_raw"]) == 4.2
