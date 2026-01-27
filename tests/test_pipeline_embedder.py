import numpy as np


def test_embedder_batch_small_model_cpu():
    from lasso_workbench.pipeline.esm_embedder import ESM2Embedder

    embedder = ESM2Embedder(model_name="facebook/esm2_t6_8M_UR50D", device="cpu")
    seqs = [
        "MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG",
        "MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG",
    ]
    emb = embedder.embed_batch(seqs)

    assert isinstance(emb, np.ndarray)
    assert emb.shape[0] == 2
    assert emb.shape[1] > 0
    assert np.isfinite(emb).all()
