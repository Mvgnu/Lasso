import numpy as np

from lasso_workbench.core.embedding_scoring import calculate_grouped_top_n_means


def test_grouped_top_n_means_locus_grouping():
    # References share a locus token; grouping should use locus= value.
    sim = np.array([[0.9, 0.2, 0.4]], dtype=np.float32)
    ref_ids = ["core|locus=L1|v1", "core|locus=L1|v2", "core|locus=L2|v1"]

    top1 = calculate_grouped_top_n_means(sim, ref_ids, top_n=1)
    assert np.allclose(top1, np.array([0.9], dtype=np.float32))

    top2 = calculate_grouped_top_n_means(sim, ref_ids, top_n=2)
    assert np.allclose(top2, np.array([0.65], dtype=np.float32), atol=1e-6)


def test_grouped_top_n_means_requires_locus():
    sim = np.array([[0.1, 0.2]], dtype=np.float32)
    ref_ids = ["core|v1", "core|v2"]
    try:
        calculate_grouped_top_n_means(sim, ref_ids, top_n=1)
    except ValueError:
        return
    raise AssertionError("Expected ValueError for missing locus token")
