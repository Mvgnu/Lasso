"""
Tests for DataService (formerly Dataset Module).

Tests CRUD operations, import/export, and file safety.
"""

import json
import tempfile
from pathlib import Path
from unittest.mock import patch
import pytest
import pandas as pd

import lasso_workbench
from lasso_workbench.services.data_service import DataService

@pytest.fixture
def temp_data_dir(tmp_path):
    """Fixture that provides a temporary data directory."""
    test_data_dir = tmp_path / "data"
    test_data_dir.mkdir()
    (test_data_dir / "training").mkdir()
    
    with patch.object(lasso_workbench, 'DATA_DIR', test_data_dir):
        with patch.object(DataService, 'DATA_DIR', test_data_dir):
            # Patch default path to use temp dir
            with patch.object(DataService, 'DEFAULT_DATASET', test_data_dir / "training" / "train.json"):
                yield test_data_dir

@pytest.fixture
def data_service(temp_data_dir):
    return DataService()


class TestLoadSaveData:
    """Tests for load/save operations."""
    
    def test_load_empty_returns_empty_df(self, data_service):
        """Test loading when no data file exists."""
        df = data_service.load_dataset()
        assert df.empty
        assert "name" in df.columns # Should have default columns
    
    def test_save_and_load_roundtrip(self, data_service):
        """Test that saved data can be loaded correctly."""
        test_data = pd.DataFrame([
            {"name": "test_1", "core": "GYPFIDACLYTGG", "source": "test", "idx": 0},
            {"name": "test_2", "core": "SLAFVPFFA", "source": "test", "idx": 1},
        ])
        
        data_service.save_dataset(test_data)
        loaded = data_service.load_dataset()
        
        assert len(loaded) == 2
        assert loaded.iloc[0]["name"] == "test_1"
        assert loaded.iloc[1]["core"] == "SLAFVPFFA"


class TestCRUDOperations:
    """Tests for add, update, delete operations."""
    
    def test_add_to_dataset(self, data_service):
        """Test adding an entry."""
        data_service.add_entry({
            "name": "Citrocin",
            "full_precursor": "MKELQTAGLSGGYPFIDACLYTGG",
            "leader": "MKELQTAGLSG",
            "core": "GYPFIDACLYTGG",
            "source": "lab"
        })
        
        df = data_service.load_dataset()
        assert len(df) == 1
        assert df.iloc[0]["name"] == "Citrocin"
    
    def test_delete_from_dataset(self, data_service):
        """Test deleting an entry."""
        data_service.add_entry({"name": "Entry1", "core": "AAAAAA", "source": "test"})
        data_service.add_entry({"name": "Entry2", "core": "BBBBBB", "source": "test"})
        
        df = data_service.load_dataset()
        assert len(df) == 2
        
        # Delete first entry (idx 0)
        success = data_service.delete_entry(0)
        assert success
        
        df = data_service.load_dataset()
        assert len(df) == 1
        assert df.iloc[0]["name"] == "Entry2"
    
    def test_update_dataset_entry(self, data_service):
        """Test updating an existing entry."""
        data_service.add_entry({"name": "Original", "core": "AAAAAA", "source": "test"})
        
        # Update it
        success = data_service.update_entry(0, {
            "name": "Updated",
            "core": "GYPFIDACLYTGG"
        })
        assert success
        
        df = data_service.load_dataset()
        assert len(df) == 1
        assert df.iloc[0]["name"] == "Updated"
        assert df.iloc[0]["core"] == "GYPFIDACLYTGG"


class TestImportOperations:
    """Tests for import operations."""

    def test_import_fasta(self, data_service, tmp_path):
        """Test importing FASTA file."""
        fasta_file = tmp_path / "test.fasta"
        fasta_file.write_text(">seq1\nMVVLAS\n>seq2\nGVVLAS")
        
        count = data_service.import_fasta(fasta_file)
        assert count == 2
        
        df = data_service.load_dataset()
        assert len(df) == 2
        assert "seq1" in df["name"].values
