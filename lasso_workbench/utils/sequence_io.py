"""
Sequence I/O helpers backed by Biopython.
"""
from __future__ import annotations

from pathlib import Path
from typing import Iterable, Iterator, List, Optional, Sequence, Tuple

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta_pairs(path: Path) -> List[Tuple[str, str]]:
    """Return (id, sequence) pairs from a FASTA file."""
    if not path.exists() or path.stat().st_size == 0:
        return []
    records: List[Tuple[str, str]] = []
    for record in SeqIO.parse(str(path), "fasta"):
        seq = str(record.seq).upper()
        if not seq:
            continue
        records.append((record.id, seq))
    return records


def write_fasta_pairs(pairs: Sequence[Tuple[str, str]], path: Path) -> None:
    """Write (id, sequence) pairs to a FASTA file."""
    records = [SeqRecord(Seq(seq), id=str(name), description="") for name, seq in pairs]
    with path.open("w") as handle:
        SeqIO.write(records, handle, "fasta")


def write_fasta_entries(entries: Iterable[Tuple[str, str, str]], path: Path) -> None:
    """Write (id, sequence, description) entries to a FASTA file."""
    records = [
        SeqRecord(Seq(seq), id=str(name), description=str(description) if description else "")
        for name, seq, description in entries
    ]
    with path.open("w") as handle:
        SeqIO.write(records, handle, "fasta")


def read_genbank_record(path: Path) -> Optional[SeqRecord]:
    """Return the first GenBank record from a file."""
    try:
        return next(SeqIO.parse(str(path), "genbank"))
    except Exception:
        return None


def iter_genbank_records(path: Path) -> Iterator[SeqRecord]:
    """Iterate GenBank records from a file."""
    return SeqIO.parse(str(path), "genbank")
