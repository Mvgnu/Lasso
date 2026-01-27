from __future__ import annotations

from pathlib import Path

from lasso_workbench.pipeline import semantic_pipeline as sp
from lasso_workbench.schemas.pipeline import RankingConfig


def test_run_semantic_pipeline_emits_stage_messages(tmp_path: Path):
    gbk_dir = Path(__file__).resolve().parents[1] / "data" / "antismash_lasso" / "gbk"
    candidates = sorted(gbk_dir.glob("*.gbk"))
    if not candidates:
        import pytest

        pytest.skip("No GBK test file found under data/antismash_lasso/gbk")
    gbk = candidates[0]

    validated = tmp_path / "validated_small.faa"
    validated.write_text(">ref1|locus=1\nMSTAVVLALAVALAGASA\n>ref2|locus=2\nMKKFFDSRREQ\n")

    messages: list[str] = []

    def cb(msg: str) -> None:
        messages.append(msg)

    # Keep it fast: small ORF window, CPU, and no need to assert scoring.
    sp.run_semantic_pipeline(
        gbk_files=[gbk],
        validated_faa=validated,
        output_dir=tmp_path,
        min_aa=80,
        max_aa=120,
        device="cpu",
        model_name="facebook/esm2_t6_8M_UR50D",
        progress=cb,
        ranking_config=RankingConfig(),
    )

    # Ensure each stage emits at least once.
    for tag in ("[1/6]", "[2/6]", "[3/6]", "[4/6]", "[5/6]", "[6/6]"):
        assert any(tag in m for m in messages), f"Missing stage message: {tag}"
