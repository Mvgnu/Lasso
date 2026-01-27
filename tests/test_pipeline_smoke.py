from pathlib import Path


def test_run_semantic_pipeline_smoke(tmp_path: Path):
    from lasso_workbench.pipeline.semantic_pipeline import run_semantic_pipeline
    from lasso_workbench.schemas.pipeline import RankingConfig

    gbk_dir = Path(__file__).resolve().parents[1] / "data" / "antismash_lasso" / "gbk"
    candidates = sorted(gbk_dir.glob("*.gbk"))
    if not candidates:
        import pytest

        pytest.skip("No GBK test file found under data/antismash_lasso/gbk")

    gbk = candidates[0]

    # Use the default precursor file
    default_precursors = (
        Path(__file__).resolve().parents[1]
        / "data"
        / "precursors"
        / "precursor_proteins_verified.faa"
    )
    if not default_precursors.exists():
        import pytest

        pytest.skip("Default precursor file not found")

    results, df = run_semantic_pipeline(
        gbk_files=[gbk],
        validated_faa=default_precursors,
        output_dir=tmp_path,
        min_aa=20,
        max_aa=60,
        device="cpu",
        model_name="facebook/esm2_t6_8M_UR50D",  # Explicit model_name
        ranking_config=RankingConfig(),
    )

    assert isinstance(results, list)
    assert df is not None
    assert "best_similarity" in df.columns
