from lasso_workbench.pipeline.semantic_pipeline import extract_orfs, build_pipeline_results
from lasso_workbench.schemas.pipeline import RankingConfig
from lasso_workbench.schemas.pipeline import BGCSegment, AnnotatedCDS, PipelineMetadata


def test_orf_coordinates_and_overlap() -> None:
    # Sequence has one ORF: ATG AAA TAA -> protein length 2, coords [0, 6)
    seq = "ATGAAATAA"
    segment = BGCSegment(
        record_id="test1",
        index=1,
        sequence=seq,
        annotated_cds=[AnnotatedCDS(start=2, end=5, strand="+")],
        metadata={"region_ranges": [{"start": 0, "end": len(seq)}]},
    )

    candidates, metadata_rows = extract_orfs([segment], min_aa=2, max_aa=3)
    assert candidates
    cand = candidates[0]
    assert cand.genomic_start == 0
    assert cand.genomic_end == 6

    similarity_rows = [
        {"candidate_id": cand.candidate_id, "best_similarity": 0.9, "best_match_id": "ref1"}
    ]
    metadata = PipelineMetadata(
        model_name="test",
        device="cpu",
        validated_faa="test",
        min_aa=2,
        max_aa=3,
        embed_batch_size=0,
        score_batch_size=0,
    )
    results = build_pipeline_results(
        similarity_rows,
        metadata_rows,
        {"test1": segment},
        ranking_config=RankingConfig(),
        pipeline_metadata=metadata,
    )
    assert results[0].candidates[0].candidate_id == cand.candidate_id


def test_orf_adjacent_cds_has_no_overlap() -> None:
    # CDS starts exactly at ORF end; should be adjacent but not overlapping.
    seq = "ATGAAATAA"
    segment = BGCSegment(
        record_id="test2",
        index=1,
        sequence=seq,
        annotated_cds=[AnnotatedCDS(start=6, end=9, strand="+")],
        metadata={"region_ranges": [{"start": 0, "end": len(seq)}]},
    )

    candidates, metadata_rows = extract_orfs([segment], min_aa=2, max_aa=3)
    cand = candidates[0]
    similarity_rows = [
        {"candidate_id": cand.candidate_id, "best_similarity": 0.9, "best_match_id": "ref1"}
    ]
    metadata = PipelineMetadata(
        model_name="test",
        device="cpu",
        validated_faa="test",
        min_aa=2,
        max_aa=3,
        embed_batch_size=0,
        score_batch_size=0,
    )
    results = build_pipeline_results(
        similarity_rows,
        metadata_rows,
        {"test2": segment},
        ranking_config=RankingConfig(),
        pipeline_metadata=metadata,
    )
    assert results[0].candidates[0].candidate_id == cand.candidate_id
