from __future__ import annotations

import hashlib
import logging
import sqlite3
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Sequence, Tuple

from Bio.Seq import Seq

from lasso_workbench.pipeline.orf_extraction import chunk_orfs
from lasso_workbench.utils.sequence_io import iter_genbank_records

logger = logging.getLogger(__name__)


@dataclass
class LassoRegion:
    record_id: str
    gbk_file: str
    region_index: int
    region_start: int
    region_end: int


def hash_peptide(peptide: str) -> str:
    return hashlib.blake2b(peptide.encode("utf-8"), digest_size=16).hexdigest()


def find_lasso_regions(record, gbk_file: str) -> List[LassoRegion]:
    regions: List[LassoRegion] = []
    record_id = record.id or gbk_file
    idx = 1
    for feat in record.features:
        if feat.type not in {"region", "cluster"}:
            continue
        products = feat.qualifiers.get("product", [])
        if not products:
            products = feat.qualifiers.get("products", [])
        if not any("lasso" in str(p).lower() for p in products):
            continue
        start = int(feat.location.start)
        end = int(feat.location.end)
        if start >= end:
            continue
        regions.append(
            LassoRegion(
                record_id=record_id,
                gbk_file=gbk_file,
                region_index=idx,
                region_start=start,
                region_end=end,
            )
        )
        idx += 1
    return regions


def init_db(conn: sqlite3.Connection, reset: bool = False) -> None:
    if reset:
        conn.executescript(
            """
            DROP TABLE IF EXISTS peptide_stats;
            DROP TABLE IF EXISTS peptide_genomes;
            DROP TABLE IF EXISTS peptide_counts;
            DROP TABLE IF EXISTS peptides;
            """
        )
    conn.executescript(
        """
        CREATE TABLE IF NOT EXISTS peptides (
            aa_hash TEXT PRIMARY KEY,
            aa_seq TEXT NOT NULL,
            aa_len INTEGER NOT NULL
        );

        CREATE TABLE IF NOT EXISTS peptide_counts (
            aa_hash TEXT PRIMARY KEY,
            instance_count INTEGER NOT NULL DEFAULT 0,
            FOREIGN KEY(aa_hash) REFERENCES peptides(aa_hash)
        );

        CREATE TABLE IF NOT EXISTS peptide_genomes (
            aa_hash TEXT NOT NULL,
            gbk_file TEXT NOT NULL,
            PRIMARY KEY (aa_hash, gbk_file),
            FOREIGN KEY(aa_hash) REFERENCES peptides(aa_hash)
        );

        CREATE INDEX IF NOT EXISTS idx_peptide_genomes_hash ON peptide_genomes(aa_hash);
        """
    )


def rebuild_stats(conn: sqlite3.Connection) -> None:
    conn.execute("DROP TABLE IF EXISTS peptide_stats")
    conn.execute(
        """
        CREATE TABLE peptide_stats AS
        SELECT
            p.aa_hash AS aa_hash,
            p.aa_len AS aa_len,
            COALESCE(g.genome_count, 0) AS genome_count,
            COALESCE(c.instance_count, 0) AS instance_count
        FROM peptides p
        LEFT JOIN (
            SELECT aa_hash, COUNT(*) AS genome_count
            FROM peptide_genomes
            GROUP BY aa_hash
        ) g ON g.aa_hash = p.aa_hash
        LEFT JOIN peptide_counts c ON c.aa_hash = p.aa_hash
        """
    )
    conn.execute("CREATE INDEX IF NOT EXISTS idx_stats_genomes ON peptide_stats(genome_count)")


def lookup_genome_counts(db_path: Path, sequences: Sequence[str], chunk_size: int = 500) -> Dict[str, int]:
    if not db_path.exists():
        return {}
    if not sequences:
        return {}

    seq_to_hash = {seq: hash_peptide(seq) for seq in sequences if seq}
    if not seq_to_hash:
        return {}

    hash_to_count: Dict[str, int] = {}
    conn = sqlite3.connect(db_path)
    try:
        has_stats = conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name='peptide_stats'"
        ).fetchone()
        if not has_stats:
            return {}
        hashes = list(set(seq_to_hash.values()))
        for i in range(0, len(hashes), chunk_size):
            batch = hashes[i:i + chunk_size]
            placeholders = ",".join(["?"] * len(batch))
            rows = conn.execute(
                f"SELECT aa_hash, genome_count FROM peptide_stats WHERE aa_hash IN ({placeholders})",
                batch,
            ).fetchall()
            for aa_hash, genome_count in rows:
                hash_to_count[aa_hash] = int(genome_count)
    finally:
        conn.close()

    return {seq: hash_to_count.get(aa_hash, 0) for seq, aa_hash in seq_to_hash.items()}


def iter_lasso_orfs(
    gbk_path: Path,
    min_aa: int,
    max_aa: int,
) -> Iterator[Tuple[LassoRegion, object]]:
    try:
        for record in iter_genbank_records(gbk_path):
            regions = find_lasso_regions(record, gbk_path.name)
            if not regions:
                continue
            seq = record.seq
            for region in regions:
                window_seq = Seq(str(seq[region.region_start:region.region_end]))
                forward_orfs = list(chunk_orfs(window_seq, "+", region.region_start, region.region_end, min_aa, max_aa))
                rc_seq = window_seq.reverse_complement()
                reverse_orfs = list(chunk_orfs(rc_seq, "-", region.region_start, region.region_end, min_aa, max_aa))
                for orf in forward_orfs + reverse_orfs:
                    yield region, orf
    except Exception as exc:  # noqa: BLE001
        logger.warning("Failed parsing %s: %s", gbk_path, exc)


def build_lasso_orf_index(
    gbk_dir: Path,
    db_path: Path,
    min_aa: int = 20,
    max_aa: int = 120,
    reset: bool = False,
    batch_size: int = 5000,
) -> dict:
    gbk_files = sorted(gbk_dir.glob("*.gbk"))
    if not gbk_files:
        raise FileNotFoundError(f"No .gbk files found in {gbk_dir}")

    conn = sqlite3.connect(db_path)
    conn.execute("PRAGMA journal_mode=WAL")
    conn.execute("PRAGMA synchronous=NORMAL")
    init_db(conn, reset=reset)

    peptide_rows: List[Tuple[str, str, int]] = []
    count_rows: List[Tuple[str, int]] = []
    genome_rows: List[Tuple[str, str]] = []

    total_orfs = 0
    for gbk_path in gbk_files:
        for region, orf in iter_lasso_orfs(gbk_path, min_aa, max_aa):
            peptide = str(orf.protein)
            if not peptide:
                continue
            aa_hash = hash_peptide(peptide)
            peptide_rows.append((aa_hash, peptide, len(peptide)))
            count_rows.append((aa_hash, 1))
            genome_rows.append((aa_hash, region.gbk_file))
            total_orfs += 1

            if len(count_rows) >= batch_size:
                _flush(conn, peptide_rows, count_rows, genome_rows)
                peptide_rows.clear()
                count_rows.clear()
                genome_rows.clear()

    if count_rows:
        _flush(conn, peptide_rows, count_rows, genome_rows)

    rebuild_stats(conn)
    conn.commit()
    conn.close()

    return {
        "gbk_files": len(gbk_files),
        "orfs_indexed": total_orfs,
        "db_path": str(db_path),
    }


def _flush(
    conn: sqlite3.Connection,
    peptide_rows: Sequence[Tuple[str, str, int]],
    count_rows: Sequence[Tuple[str, int]],
    genome_rows: Sequence[Tuple[str, str]],
) -> None:
    conn.executemany(
        "INSERT OR IGNORE INTO peptides (aa_hash, aa_seq, aa_len) VALUES (?, ?, ?)",
        peptide_rows,
    )
    conn.executemany(
        """
        INSERT INTO peptide_counts (aa_hash, instance_count)
        VALUES (?, ?)
        ON CONFLICT(aa_hash) DO UPDATE SET instance_count = instance_count + excluded.instance_count
        """,
        count_rows,
    )
    conn.executemany(
        "INSERT OR IGNORE INTO peptide_genomes (aa_hash, gbk_file) VALUES (?, ?)",
        genome_rows,
    )
