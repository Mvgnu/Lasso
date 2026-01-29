import random

from lasso_workbench.altframe.codon_shuffle import shuffle_synonymous_codons, synonym_maps


def test_shuffle_preserves_start_and_multiset():
    rng = random.Random(42)
    codons = ["ATG", "GCT", "GCC", "GCA", "GCG"]
    aas = ["M", "A", "A", "A", "A"]
    shuffled = shuffle_synonymous_codons(codons, aas, rng)

    assert shuffled[0] == "ATG"
    assert sorted(shuffled[1:]) == sorted(codons[1:])

    _, codon_to_aa = synonym_maps()
    assert [codon_to_aa.get(c, "X") for c in shuffled] == aas


def test_shuffle_does_not_move_unknown_aa():
    rng = random.Random(7)
    codons = ["ATG", "GCT", "NNN", "GCC"]
    aas = ["M", "A", "X", "A"]
    shuffled = shuffle_synonymous_codons(codons, aas, rng)

    assert shuffled[0] == "ATG"
    assert shuffled[2] == "NNN"


def test_shuffle_deterministic_with_seed():
    codons = ["ATG", "GCT", "GCC", "GCA", "GCG"]
    aas = ["M", "A", "A", "A", "A"]
    rng1 = random.Random(13)
    rng2 = random.Random(13)
    assert shuffle_synonymous_codons(codons, aas, rng1) == shuffle_synonymous_codons(codons, aas, rng2)


def test_shuffle_preserves_first_position_even_with_synonyms():
    rng = random.Random(5)
    codons = ["GTG", "GTT", "GTC", "GTA"]
    aas = ["V", "V", "V", "V"]
    shuffled = shuffle_synonymous_codons(codons, aas, rng)
    assert shuffled[0] == "GTG"


def test_shuffle_preserves_per_aa_codons():
    rng = random.Random(9)
    codons = ["ATG", "CTG", "TTA", "TCT", "AGC", "CGT", "AGA"]
    aas = ["M", "L", "L", "S", "S", "R", "R"]
    shuffled = shuffle_synonymous_codons(codons, aas, rng)

    per_aa_before = {
        "L": sorted([codons[1], codons[2]]),
        "S": sorted([codons[3], codons[4]]),
        "R": sorted([codons[5], codons[6]]),
    }
    per_aa_after = {
        "L": sorted([shuffled[1], shuffled[2]]),
        "S": sorted([shuffled[3], shuffled[4]]),
        "R": sorted([shuffled[5], shuffled[6]]),
    }
    assert per_aa_before == per_aa_after
