from lasso_workbench.core.prediction import RuleEngine, load_active_ruleset


def test_score_cleavage_site_slicing_boundaries():
    rules = load_active_ruleset()
    engine = RuleEngine(rules)
    sequence = "MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG"

    # Position 0 -> empty leader, full core.
    score0, _ = engine.score_cleavage_site(sequence, 0)
    assert isinstance(score0, float)

    # Position at end -> full leader, empty core.
    score_end, _ = engine.score_cleavage_site(sequence, len(sequence))
    assert isinstance(score_end, float)

    # A mid-position should split cleanly.
    pos = 10
    leader = sequence[:pos]
    core = sequence[pos:]
    assert leader + core == sequence
