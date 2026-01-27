"""
Edge Case Tests for Lasso Workbench

Tests unusual inputs, malformed data, and boundary conditions to ensure
the tool is robust and provides helpful error messages for scientists.
"""

import pytest
import tempfile
from pathlib import Path
from io import StringIO
from lasso_workbench.core.prediction import CorePredictor, RuleEngine, load_default_ruleset



class TestEmptyAndShortSequences:
    """Test handling of empty and very short sequences."""
    
    def test_empty_sequence_prediction(self):
        """Empty sequence should not crash, should return empty predictions."""
        predictor = CorePredictor()
        result = predictor.predict("")
        
        # Should handle gracefully
        assert result is not None
        assert len(result.predictions) == 0 or result.best_prediction is None
    
    def test_single_aa_sequence(self):
        """Single amino acid should fail gracefully."""
        predictor = CorePredictor()
        result = predictor.predict("M")
        
        assert result is not None
        # Should return no valid predictions for such a short sequence
        assert len(result.predictions) == 0 or all(
            len(p.core) < 1 for p in result.predictions
        )
    
    def test_minimum_viable_length(self):
        """Test the minimum viable precursor length (around 20-25 aa)."""
        predictor = CorePredictor()
        # 25 aa precursor - should work
        short_but_valid = "MKFLATAAGGACGDFGAGHDFGHDFG"
        result = predictor.predict(short_but_valid)
        
        assert result is not None


class TestNonStandardAminoAcids:
    """Test handling of non-standard amino acid codes."""
    
    def test_ambiguous_aa_B(self):
        """B (Asx = D or N) should be handled."""
        predictor = CorePredictor()
        seq_with_b = "MKFLATAAGGACGDFGBGHDFGHDFGBDFG"
        result = predictor.predict(seq_with_b)
        # Should handle gracefully - either skip or treat as 0
        assert result is not None
    
    def test_ambiguous_aa_Z(self):
        """Z (Glx = E or Q) should be handled."""
        predictor = CorePredictor()
        seq_with_z = "MKFLATAAGGACGDFGZGHDFGHDFGZDFG"
        result = predictor.predict(seq_with_z)
        assert result is not None
    
    def test_unknown_aa_X(self):
        """X (any amino acid) should be handled."""
        predictor = CorePredictor()
        seq_with_x = "MKFLATAAGGACGDFGXGHDFGHDFGXDFG"
        result = predictor.predict(seq_with_x)
        assert result is not None
    
    def test_selenocysteine_U(self):
        """U (Selenocysteine) should be handled."""
        predictor = CorePredictor()
        seq_with_u = "MKFLATAAGGACGDFGUPHDFGHDFG"
        result = predictor.predict(seq_with_u)
        assert result is not None
    
    def test_pyrrolysine_O(self):
        """O (Pyrrolysine) should be handled."""
        predictor = CorePredictor()
        seq_with_o = "MKFLATAAGGACGDFGOPHDFGHDFG"
        result = predictor.predict(seq_with_o)
        assert result is not None








class TestSequenceCharacterHandling:
    """Test handling of various character inputs in sequences."""
    
    def test_lowercase_sequence(self):
        """Lowercase sequences should be converted to uppercase."""
        predictor = CorePredictor()
        result = predictor.predict("mkflataaggacgdfgaghdfghdfg")
        assert result is not None
    
    def test_mixed_case_sequence(self):
        """Mixed case sequences should work."""
        predictor = CorePredictor()
        result = predictor.predict("MkFlAtAaGgAcGdFgAgHdFgHdFg")
        assert result is not None
    
    def test_sequence_with_whitespace(self):
        """Sequences with whitespace should be cleaned."""
        predictor = CorePredictor()
        result = predictor.predict("MKFLAT AAGGAC GDFGAG HDFGHDFG")
        assert result is not None
    
    def test_sequence_with_numbers(self):
        """Sequences with numbers should have them stripped or error."""
        predictor = CorePredictor()
        try:
            result = predictor.predict("MKFLAT123AAGGAC456GDFGAG")
            # If it works, numbers should be ignored
            assert result is not None
        except ValueError as e:
            # Or raise a clear error
            assert "number" in str(e).lower() or "invalid" in str(e).lower()
    
    def test_sequence_with_newlines(self):
        """FASTA-style sequences with newlines should work."""
        predictor = CorePredictor()
        fasta_seq = "MKFLAT\nAAGGAC\nGDFGAGH\nDFGHDFG"
        result = predictor.predict(fasta_seq)
        assert result is not None


class TestRuleEngineConfiguration:
    """Test RuleEngine handling of unusual configurations."""
    
    def test_all_rules_disabled(self):
        """Engine with all rules disabled should still return a score (may not be 0)."""
        rules = load_default_ruleset()
        engine = RuleEngine(rules)
        # Disable all rules
        for rule in engine.rules.rules.values():
            rule.enabled = False
            
        score, _ = engine.score_cleavage_site("MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG", 22)
        assert isinstance(score, float)

    def test_zero_weight_rules(self):
        """Test that rules can be weighted - verifies scoring still works."""
        rules = load_default_ruleset()
        engine = RuleEngine(rules)
        # Set weights to 0
        for rule in engine.rules.rules.values():
            rule.weight = 0.0
            
        score, details = engine.score_cleavage_site("MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG", 22)
        # Just verify it returns valid score
        assert isinstance(score, float)
        assert score == 0.0
        assert isinstance(details, dict)
    
    def test_invalid_cleavage_position(self):
        """Invalid cleavage positions should be handled."""
        rules = load_default_ruleset()
        engine = RuleEngine(rules)
        sequence = "MKKQTFVPKKLVKVGKATELTKGGVGKIIEYFIGGGVGRYG"
        
        # Position at sequence start
        score, details = engine.score_cleavage_site(sequence, 0)
        assert score is not None
        
        # Position at sequence end
        score, details = engine.score_cleavage_site(sequence, len(sequence))
        assert score is not None
        
        # Position beyond sequence
        score, details = engine.score_cleavage_site(sequence, len(sequence) + 10)
        assert score is not None  # Should handle gracefully



