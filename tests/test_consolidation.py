"""Tests for pipeline consolidation helpers."""

import numpy as np
import pytest


class TestCosineSimilarity:
    """Test the sklearn-based cosine similarity helper."""

    def test_identical_vectors_have_similarity_one(self):
        """Identical vectors should have similarity 1.0."""
        from sklearn.metrics.pairwise import cosine_similarity
        
        vec = np.array([[1.0, 2.0, 3.0]])
        sims = cosine_similarity(vec, vec, dense_output=True)
        
        assert sims.shape == (1, 1)
        np.testing.assert_almost_equal(sims[0, 0], 1.0, decimal=5)

    def test_orthogonal_vectors_have_similarity_zero(self):
        """Orthogonal vectors should have similarity 0.0."""
        from sklearn.metrics.pairwise import cosine_similarity
        
        vec1 = np.array([[1.0, 0.0, 0.0]])
        vec2 = np.array([[0.0, 1.0, 0.0]])
        sims = cosine_similarity(vec1, vec2, dense_output=True)
        
        assert sims.shape == (1, 1)
        np.testing.assert_almost_equal(sims[0, 0], 0.0, decimal=5)

    def test_opposite_vectors_have_similarity_minus_one(self):
        """Opposite vectors should have similarity -1.0."""
        from sklearn.metrics.pairwise import cosine_similarity
        
        vec1 = np.array([[1.0, 2.0, 3.0]])
        vec2 = np.array([[-1.0, -2.0, -3.0]])
        sims = cosine_similarity(vec1, vec2, dense_output=True)
        
        assert sims.shape == (1, 1)
        np.testing.assert_almost_equal(sims[0, 0], -1.0, decimal=5)

    def test_batch_similarity(self):
        """Test batch similarity computation."""
        from sklearn.metrics.pairwise import cosine_similarity
        
        candidates = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [1.0, 1.0, 0.0],
        ])
        references = np.array([
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
        ])
        
        sims = cosine_similarity(candidates, references, dense_output=True)
        
        assert sims.shape == (3, 2)
        # First candidate is identical to first reference
        np.testing.assert_almost_equal(sims[0, 0], 1.0, decimal=5)
        # First candidate is orthogonal to second reference
        np.testing.assert_almost_equal(sims[0, 1], 0.0, decimal=5)


class TestTranslateBacterial:
    """Test bacterial translation with alternative start codons."""

    def test_gtg_start_becomes_methionine(self):
        """GTG start codon should translate to M, not V."""
        from lasso_workbench.utils.translation import translate_bacterial
        
        dna = "GTGAAATAA"  # GTG-AAA-TAA = (V|M)-K-*
        protein = translate_bacterial(dna)
        
        assert protein.startswith("M"), f"Expected M, got {protein[0]}"
        assert protein == "MK"

    def test_ttg_start_becomes_methionine(self):
        """TTG start codon should translate to M, not L."""
        from lasso_workbench.utils.translation import translate_bacterial
        
        dna = "TTGAAATAG"  # TTG-AAA-TAG = (L|M)-K-*
        protein = translate_bacterial(dna)
        
        assert protein.startswith("M"), f"Expected M, got {protein[0]}"
        assert protein == "MK"

    def test_atg_start_unchanged(self):
        """ATG start codon should stay as M."""
        from lasso_workbench.utils.translation import translate_bacterial
        
        dna = "ATGAAATAA"  # ATG-AAA-TAA = M-K-*
        protein = translate_bacterial(dna)
        
        assert protein == "MK"

    def test_standard_translation_would_be_wrong(self):
        """Verify standard BioPython translate gives wrong result for GTG."""
        from Bio.Seq import Seq
        
        dna = "GTGAAATAA"
        standard = str(Seq(dna).translate(table=11, to_stop=True))
        
        # BioPython gives V, which is wrong for start codon
        assert standard.startswith("V"), "BioPython should give V (which is wrong)"

    def test_internal_gtg_remains_valine(self):
        """GTG within the sequence (not at start) should remain Valine."""
        from lasso_workbench.utils.translation import translate_bacterial

        # ATG-GTG-TAA = M-V-*
        dna = "ATGGTGTAA"
        protein = translate_bacterial(dna)
        
        assert protein == "MV"



class TestBioPythonConsolidation:
    """Test that BioPython is used for parsing functions."""

    def test_fasta_import_uses_biopython(self):
        """FASTA reading should use BioPython SeqIO."""
        import inspect
        from lasso_workbench.utils.sequence_io import read_fasta_pairs
        
        source = inspect.getsource(read_fasta_pairs)
        assert "SeqIO.parse" in source

    def test_genbank_read_uses_biopython(self):
        """GenBank reading should use BioPython SeqIO."""
        import inspect
        from lasso_workbench.utils.sequence_io import read_genbank_record
        
        source = inspect.getsource(read_genbank_record)
        assert "SeqIO.parse" in source

