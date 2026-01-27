"""Tests for service layer."""

import pytest
import shutil
from pathlib import Path
import tempfile
import json
import pandas as pd

from lasso_workbench.services.rules_service import RulesService
from lasso_workbench.services.data_service import DataService
from lasso_workbench.schemas.config import RuleSet

class TestRulesService:
    def test_load_default_rules(self):
        service = RulesService()
        rules = service.get_rules()
        assert isinstance(rules, RuleSet)
        assert "core_start_residue" in rules.rules

    def test_save_and_load_user_rules(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "user_rules.json"
            service = RulesService(config_path=path)
            
            # Modify a rule
            rules = service.get_rules()
            if "core_start_residue" in rules.rules:
                 service.update_rule("core_start_residue", {"weight": 99.9})
            
            # Save
            service.save_rules()
            
            # Reload
            service2 = RulesService(config_path=path)
            rules2 = service2.get_rules()
            assert rules2.rules["core_start_residue"].weight == 99.9


class TestDataService:
    def test_init_creates_dir(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "data" / "dataset.csv"
            service = DataService(data_path=path)
            assert path.parent.exists()

    def test_add_and_list_entry(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "dataset.csv"
            service = DataService(data_path=path)
            
            entry = {
                "name": "Test Entry",
                "sequence": "MKKQTF",
                "core": "TF",
                "source": "test"
            }
            service.add_entry(entry)
            
            df = service.load_dataset()
            assert len(df) == 1
            assert df.iloc[0]["name"] == "Test Entry"

    def test_update_entry(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "dataset.csv"
            service = DataService(data_path=path)
            
            service.add_entry({"name": "Old Name", "sequence": "ABC"})
            service.update_entry(0, {"name": "New Name"})
            
            df2 = service.load_dataset()
            assert df2.iloc[0]["name"] == "New Name"

    def test_delete_entry(self):
        with tempfile.TemporaryDirectory() as td:
            path = Path(td) / "dataset.csv"
            service = DataService(data_path=path)
            
            service.add_entry({"name": "To Delete"})
            service.delete_entry(0)
            
            df2 = service.load_dataset()
            assert len(df2) == 0

