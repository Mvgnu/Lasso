from Bio.Seq import Seq
from lasso_workbench.pipeline.orf_extraction import chunk_orfs
from lasso_workbench.schemas.pipeline import CandidateORF, BGCSegment, AnnotatedCDS, TransientORF
from lasso_workbench.pipeline.semantic_pipeline import extract_orfs

def test_chunk_orfs_multi_start():
    # Sequence: ATG(M)...GTG(M)...TAG(*)
    # Counts: 
    # Frame 0: ATG...TAG (Length depends on sequence)
    # If we have ATG (start 0) ... GTG (start 9) ... TAG (stop 18)
    # Sequence: ATGGGGGGGGTGGGGGGGTAG
    # Indices:  012345678901234567890
    # Pos 0: ATG
    # Pos 9: GTG (valid start in bacterial)
    # Pos 18: TAG (stop)
    
    seq = Seq("ATGGGGGGGGTGGGGGGGTAG") 
    # Length: 21 nt.
    # ORF 1: 0-18 (18nt, 6aa)
    # ORF 2: 9-18 (9nt, 3aa)
    
    # min_aa=3 to catch both
    orfs = list(chunk_orfs(seq, "+", 1, 21, min_aa=3, max_aa=100))
    
    # Should get 2 ORFs
    assert len(orfs) == 2
    
    # Verify starts
    starts = sorted([orf.start_codon for orf in orfs])
    assert starts == ["ATG", "GTG"]

def test_extract_orfs_basic():
    seq = Seq("ATGGGGGGGGTGGGGGGGTAG") 
    segment = BGCSegment(
        record_id="test_rec",
        index=1,
        sequence=str(seq),
        annotated_cds=[],
        metadata={"region_ranges": [{"start": 0, "end": len(seq)}]},
    )
    
    candidates, _ = extract_orfs([segment], min_aa=3, max_aa=100)
    assert len(candidates) >= 2
    
    # Check candidate IDs are unique
    ids = {c.candidate_id for c in candidates}
    assert len(ids) == len(candidates)

def test_schema_candidate_orf_basic():
    # Basic schema construction should still work.
    cand = CandidateORF(
        candidate_id="id",
        record_id="rec",
        strand="+",
        frame=0,
        aa_length=10,
        nt_length=30,
        protein_sequence="MMMMMMMMMM",
        dna_sequence="ATG" * 10,
        genomic_start=1,
        genomic_end=30,
    )
    assert cand.candidate_id == "id"
