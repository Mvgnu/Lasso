"""
ORF extraction utilities for 6-frame translation.

Extracts all potential open reading frames (ORFs) from a DNA sequence,
handling bacterial translation rules (alternative start codons).

Instead of using the existing ORF finders which filter on assumptions, 
we implement our own version
"""
from __future__ import annotations

from typing import Iterator

from Bio.Seq import Seq
from lasso_workbench.schemas.pipeline import TransientORF
from lasso_workbench.utils.translation import (
    translate_bacterial,
    START_CODONS,
    BACTERIAL_CODON_TABLE,
)

def chunk_orfs(
    sequence: Seq,
    strand: str,
    window_start: int,
    window_end: int,
    min_aa: int,
    max_aa: int
) -> Iterator[TransientORF]:
    """
    Extract ORFs from a sequence window on one strand.

    Strategy:
    - For each frame, translate the full frame (bacterial table) to find stop codons.
    - Collect all start codons in that frame.
    - For each stop, emit ORFs for all starts since the previous stop
      (keeps overlapping ORFs by design).
    
    Args:
        sequence: Bio.Seq object (nucleotide sequence)
        strand: Strand identifier ('+' or '-')
        window_start: Genomic start coordinate of the window
        window_end: Genomic end coordinate of the window
        min_aa: Minimum amino acid length
        max_aa: Maximum amino acid length
        
    Yields:
        TransientORF objects
    """
    seq_str = str(sequence).upper()
    for frame in (0, 1, 2):
        # Frame-specific DNA: shift by 0/1/2 nt so codons align at indices 0,3,6,...
        frame_seq = seq_str[frame:]
        if len(frame_seq) < 3:
            continue
        
        # Trim to a length divisible by 3 so translation does not include a partial codon.
        frame_len = len(frame_seq) - (len(frame_seq) % 3)
        frame_seq = frame_seq[:frame_len]
        if frame_len < 3:
            continue

        # Get all start codon positions (nt offsets within this frame_seq).
        # These offsets are multiples of 3.
        start_positions = [
            offset
            for offset in range(0, frame_len, 3)
            if frame_seq[offset:offset + 3] in START_CODONS
        ]
        if not start_positions:
            continue

        # Translate full frame to find stop codons ("*").
        protein = str(Seq(frame_seq).translate(table=BACTERIAL_CODON_TABLE))
        # Translate the entire frame in one go to locate stops.
        # '*' marks stop codons. Position in protein -> map back to DNA via *3.
        stop_positions = [idx * 3 for idx, aa in enumerate(protein) if aa == "*"]
        
        if not stop_positions:
            continue

        # start_idx tracks the first start >= current stop; everything before it is < stop.
        start_idx = 0 

        for stop_nt in stop_positions:
            # Keep track of the previous start_idx to avoid re-processing the same starts.
            prev_start_idx = start_idx
            
            # Advance start_idx to include all starts strictly before this stop by comparing 
            # start_positions[start_idx] (nt offset of start codon) with stop_nt (nt offset of stop codon)
            # Those starts define ORFs that terminate at this stop codon.
            while start_idx < len(start_positions) and start_positions[start_idx] < stop_nt:
                start_idx += 1

            # No new starts before this stop: this stop yields no ORFs.
            if start_idx == prev_start_idx:
                continue

            # The stop codon position in the original sequence (accounting for frame offset).
            nt_end = frame + stop_nt

            # Emit one ORF for each start codon that ends at this stop. 
            for start_offset in start_positions[prev_start_idx:start_idx]:
                nt_start = frame + start_offset
                nt_len = nt_end - nt_start
                aa_len = nt_len // 3

                if not (min_aa <= aa_len <= max_aa):
                    continue

                dna_slice = sequence[nt_start:nt_end]
                if len(dna_slice) < 3:
                    continue

                # Bacterial translation: GTG/TTG → M at first position.
                protein_seq = translate_bacterial(str(dna_slice), start_codon_to_met=True)

                # Convert local window coords to genomic coords.
                if strand == "+":
                    genomic_start = window_start + nt_start
                    genomic_end = window_start + nt_end
                    frame_id = frame
                else:
                    # Reverse complement: flip coordinates. [L - end, L - start)
                    genomic_start = window_end - nt_end
                    genomic_end = window_end - nt_start
                    frame_id = -(frame + 1)
                
                yield TransientORF(
                    dna=str(dna_slice),
                    protein=protein_seq,
                    strand=strand,
                    frame=frame_id,
                    genomic_start=genomic_start,
                    genomic_end=genomic_end,
                    aa_len=aa_len,
                    nt_len=nt_len,
                    start_codon=frame_seq[start_offset:start_offset + 3],
                )

            # start_idx already points past all starts we just processed.
