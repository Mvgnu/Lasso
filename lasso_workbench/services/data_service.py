"""
DataService - file-based dataset management.
Simply: dataset file → consumers. No timestamps, no history.
"""
import logging
from pathlib import Path
from typing import Dict, Any, Optional

import pandas as pd

from lasso_workbench.utils.sequence_io import read_fasta_pairs
logger = logging.getLogger(__name__)


class DataService:
    """
    Simple file-based dataset service.
    Load, save, add, update, delete. That's it.
    """
    DATA_DIR = Path(__file__).resolve().parents[2] / "data"
    DEFAULT_DATASET = DATA_DIR / "training" / "train.json"
    PRECURSORS_DIR = DATA_DIR / "precursors"
    CURATED_PRECURSOR_TSV = PRECURSORS_DIR / "precursor_proteins_curated.tsv"

    DATASET_EXTENSIONS = {".csv", ".tsv", ".json"}
    FASTA_EXTENSIONS = {".fasta", ".fa", ".faa"}

    COLUMNS = [
        "name",
        "sequence",
        "full_precursor",
        "leader",
        "core",
        "source",
        "validation_status",
    ]
    
    def __init__(self, data_path: Optional[Path] = None):
        self.data_path = data_path or self.DEFAULT_DATASET
        self.data_path.parent.mkdir(parents=True, exist_ok=True)
        self._read_only = self._is_read_only_path(self.data_path)

    def set_dataset_path(self, path: Path) -> None:
        """Point the service at a new dataset file."""
        self.data_path = path
        self.data_path.parent.mkdir(parents=True, exist_ok=True)
        self._read_only = self._is_read_only_path(self.data_path)

    def is_read_only(self) -> bool:
        return self._read_only

    def _is_read_only_path(self, path: Path) -> bool:
        try:
            resolved = path.resolve()
        except Exception:
            resolved = path
        if self.PRECURSORS_DIR in resolved.parents:
            return resolved != self.CURATED_PRECURSOR_TSV
        return False

    def _assert_writable(self) -> None:
        if self._read_only:
            raise PermissionError(f"Dataset is read-only: {self.data_path}")

    def _read_dataset(self, path: Path) -> pd.DataFrame:
        suffix = path.suffix.lower()
        if suffix == ".json":
            return pd.read_json(path)
        if suffix == ".tsv":
            return pd.read_csv(path, sep="\t", on_bad_lines="skip")
        return pd.read_csv(path, on_bad_lines="skip")

    def _write_dataset(self, df: pd.DataFrame, path: Path) -> None:
        self._assert_writable()
        suffix = path.suffix.lower()
        if suffix == ".json":
            df.to_json(path, orient="records", indent=2)
            return
        if suffix == ".tsv":
            df.to_csv(path, sep="\t", index=False)
            return
        df.to_csv(path, index=False)

    def _merge_entries(self, new_df: pd.DataFrame, dedupe_sequence: bool = False) -> int:
        if new_df.empty:
            return 0
        df = self.load_dataset()
        combined = pd.concat([df, new_df], ignore_index=True)
        if dedupe_sequence and "sequence" in combined.columns:
            combined = combined.drop_duplicates(subset=["sequence"], keep="first")
        added = max(0, len(combined) - len(df))
        self.save_dataset(combined)
        return added

    def load_dataset(self) -> pd.DataFrame:
        """Load dataset file as DataFrame."""
        if not self.data_path.exists():
            return pd.DataFrame(columns=self.COLUMNS)
        
        try:
            df = self._read_dataset(self.data_path)
            return df
        except Exception as e:
            logger.warning("Error loading %s: %s", self.data_path, e)
            return pd.DataFrame(columns=self.COLUMNS)

    def save_dataset(self, df: pd.DataFrame) -> None:
        """Save DataFrame to dataset file."""
        self._write_dataset(df, self.data_path)

    def add_entry(self, entry: Dict[str, Any]) -> None:
        """Add a single entry."""
        new_row = pd.DataFrame([entry])
        self._merge_entries(new_row, dedupe_sequence=False)

    def update_entry(self, idx: int, updates: Dict[str, Any]) -> bool:
        """Update an entry by index."""
        df = self.load_dataset()
        if idx not in df.index:
            return False
        
        update_keys = [key for key in updates if key in df.columns]
        if update_keys:
            df.loc[idx, update_keys] = [updates[key] for key in update_keys]
        self.save_dataset(df)
        return True

    def delete_entry(self, idx: int) -> bool:
        """Delete an entry by index."""
        df = self.load_dataset()
        if idx not in df.index:
            return False
        
        df = df.drop(index=idx).reset_index(drop=True)
        self.save_dataset(df)
        return True

    def import_csv(self, file_path: Path) -> int:
        """Import entries from another CSV. Returns count added."""
        if self._read_only:
            logger.warning("Refusing import into read-only dataset: %s", self.data_path)
            return 0
        try:
            imported = pd.read_csv(file_path, on_bad_lines="skip")
            if imported.empty:
                return 0
            return self._merge_entries(imported, dedupe_sequence=True)
        except Exception as e:
            logger.warning("Failed to import %s: %s", file_path, e)
            return 0

    def load_as_list(self) -> list:
        """Load dataset as list of dicts (for rule evaluation)."""
        df = self.load_dataset()
        return df.to_dict(orient="records")

    def import_fasta(self, file_path: Path) -> int:
        """Import sequences from FASTA. Returns count added."""
        if self._read_only:
            logger.warning("Refusing import into read-only dataset: %s", self.data_path)
            return 0
        try:
            records = read_fasta_pairs(file_path)
            if not records:
                return 0
            entries = [
                {"name": record_id, "sequence": sequence, "source": "fasta"}
                for record_id, sequence in records
            ]
            return self._merge_entries(pd.DataFrame(entries), dedupe_sequence=True)
        except Exception as e:
            logger.warning("Failed to import FASTA %s: %s", file_path, e)
        return 0

    def import_file(self, file_path: Path) -> int:
        """Import entries from CSV/TSV/JSON/FASTA. Returns count added."""
        if self._read_only:
            logger.warning("Refusing import into read-only dataset: %s", self.data_path)
            return 0
        suffix = file_path.suffix.lower()
        if suffix in self.FASTA_EXTENSIONS:
            return self.import_fasta(file_path)
        if suffix in self.DATASET_EXTENSIONS:
            try:
                imported = self._read_dataset(file_path)
                if imported.empty:
                    return 0
                return self._merge_entries(imported, dedupe_sequence=True)
            except Exception as e:
                logger.warning("Failed to import %s: %s", file_path, e)
                return 0
        raise ValueError(f"Unsupported file type: {file_path.suffix}")
