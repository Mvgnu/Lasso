from lasso_workbench.core.prediction import RuleEngine, load_default_ruleset


def _score(sequence: str, pos: int):
    engine = RuleEngine(load_default_ruleset())
    total, breakdown = engine.score_cleavage_site(sequence, pos)
    return total, breakdown


def test_rule_core_start_residue_positive():
    leader = "AAAAAAAAAA"
    core = "G" + "A" * 13
    _, breakdown = _score(leader + core, len(leader))
    assert "core_start" in breakdown
    assert round(breakdown["core_start"], 2) == 1.2


def test_rule_ring_closure_residue_positive():
    leader = "AAAAAAAAAA"
    core = "AAAAAA" + "D" + "AAAAAAA"
    _, breakdown = _score(leader + core, len(leader))
    assert "ring_closure" in breakdown
    assert round(breakdown["ring_closure"], 2) == 4.5


def test_rule_steric_lock_region_positive():
    leader = "AAAAAAAAAA"
    core = "AAAAAA" + "D" + "AAA" + "F" + "A" * 9  # length 20
    _, breakdown = _score(leader + core, len(leader))
    assert "steric_lock" in breakdown
    assert round(breakdown["steric_lock"], 2) == 1.0


def test_rule_leader_penultimate_residue_positive():
    leader = "AAAAAAAATA"
    core = "A" * 14
    _, breakdown = _score(leader + core, len(leader))
    assert "leader_penultimate" in breakdown
    assert round(breakdown["leader_penultimate"], 2) == 1.0


def test_rule_leader_signature_motif_positive():
    leader = "AAAAAYQQPALAAAA"
    core = "A" * 14
    _, breakdown = _score(leader + core, len(leader))
    assert "leader_signature" in breakdown
    assert round(breakdown["leader_signature"], 2) == 1.2


def test_rule_cleavage_motif_positive():
    leader = "AAAAAAAAAG"
    core = "G" + "A" * 13
    _, breakdown = _score(leader + core, len(leader))
    assert "cleavage_motif" in breakdown
    assert round(breakdown["cleavage_motif"], 2) == 0.48


def test_rule_core_length_optimal():
    leader = "AAAAAAAAAA"
    core = "A" * 20
    _, breakdown = _score(leader + core, len(leader))
    assert "core_length" in breakdown
    assert round(breakdown["core_length"], 2) == 1.0


def test_rule_hydrophobicity_transition_positive():
    leader = "EEEEEEEEEE"
    core = "VVVVVVVVVV"
    _, breakdown = _score(leader + core, len(leader))
    assert "hydrophobicity" in breakdown
    assert breakdown["hydrophobicity"] > 0


def test_rule_charge_transition_positive():
    leader = "KKKKKKKKKK"
    core = "AAAAAAAAAA"
    _, breakdown = _score(leader + core, len(leader))
    assert "charge_transition" in breakdown
    assert breakdown["charge_transition"] > 0
