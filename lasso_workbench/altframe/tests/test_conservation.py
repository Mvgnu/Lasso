import math

from lasso_workbench.altframe.conservation import (
    mean_pairwise_identity,
    p_value,
    summarize_null,
    translate_window,
)


def test_mean_pairwise_identity_identical():
    assert mean_pairwise_identity(["MAG", "MAG", "MAG"]) == 1.0


def test_mean_pairwise_identity_different():
    result = mean_pairwise_identity(["MAG", "MAR"])
    assert abs(result - (2 / 3)) < 0.01


def test_mean_pairwise_identity_empty():
    assert math.isnan(mean_pairwise_identity([]))


def test_mean_pairwise_identity_single_sequence():
    assert math.isnan(mean_pairwise_identity(["MAG"]))


def test_mean_pairwise_identity_length_mismatch():
    assert mean_pairwise_identity(["AAAA", "AAAAZZ"]) == 1.0


def test_translate_window_frames():
    base = "ATGAAATTTCCC"
    assert translate_window(base, "+", 0, "zero").startswith("MKF")
    assert translate_window("A" + base, "+", 1, "zero").startswith("MKF")
    assert translate_window("AA" + base, "+", 2, "zero").startswith("MKF")


def test_translate_window_reverse_strand():
    seq = "ATGAAATTTCCC"
    pep = translate_window(seq, "-", -1, "zero")
    assert pep == "GKFH"


def test_translate_window_reverse_strand_frame_offset():
    seq = "ATGAAATTTCCC"
    pep = translate_window(seq, "-", -2, "zero")
    assert pep == "GNF"


def test_translate_window_terminal_stop_removed():
    seq = "ATGAAATAA"
    pep = translate_window(seq, "+", 0, "zero")
    assert pep == "MK"


def test_translate_window_internal_stop_preserved():
    seq = "ATGTAAATG"
    pep = translate_window(seq, "+", 0, "zero")
    assert "*" in pep

def test_p_value_basic():
    assert p_value(1, 3) == 0.5


def test_summarize_null_basic():
    mean, std = summarize_null([0.2, 0.4, 0.6])
    assert abs(mean - 0.4) < 1e-6
    assert std > 0
