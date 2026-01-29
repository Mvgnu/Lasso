from lasso_workbench.altframe.binning import compute_bins, window_from_bins


def test_compute_bins_basic():
    result = compute_bins(1000, 2000, 1150, 1600, 20)
    assert result == (3, 11)


def test_compute_bins_no_overlap():
    result = compute_bins(1000, 2000, 500, 800, 20)
    assert result is None


def test_compute_bins_clamps_end():
    result = compute_bins(0, 1000, 900, 1000, 20)
    assert result is not None
    bin_start, bin_end = result
    assert bin_end == 19
    assert bin_start <= bin_end


def test_compute_bins_one_nt_overlap():
    result = compute_bins(0, 1000, 999, 1000, 20)
    assert result == (19, 19)


def test_window_from_bins_roundtrip_includes_original():
    gene_start = 0
    gene_end = 1000
    bin_start = 4
    bin_end = 6
    bins = 20
    window_start, window_end = window_from_bins(gene_start, gene_end, bin_start, bin_end, bins)
    roundtrip = compute_bins(gene_start, gene_end, window_start, window_end, bins)
    assert roundtrip == (bin_start, bin_end)


def test_window_roundtrip_randomized():
    import random

    rng = random.Random(0)
    for _ in range(50):
        gene_start = 0
        gene_end = rng.randint(50, 2000)
        bins = rng.randint(5, 50)
        bin_start = rng.randint(0, bins - 1)
        bin_end = rng.randint(bin_start, bins - 1)
        window_start, window_end = window_from_bins(gene_start, gene_end, bin_start, bin_end, bins)
        roundtrip = compute_bins(gene_start, gene_end, window_start, window_end, bins)
        assert roundtrip == (bin_start, bin_end)
