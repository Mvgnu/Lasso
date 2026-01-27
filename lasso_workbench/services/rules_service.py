import logging
from pathlib import Path
from typing import Optional, Dict, Any

from lasso_workbench.schemas.config import RuleSet, RuleParameter
from lasso_workbench.config.rules_io import (
    DEFAULT_RULES_PATH,
    USER_RULES_PATH,
    PRESET_DIR,
    load_ruleset,
    load_active_ruleset,
    save_ruleset,
    list_presets,
)

logger = logging.getLogger(__name__)

class RulesService:
    """
    Service for managing rule configurations.
    Handles loading, saving, and updating rule sets.
    """
    DEFAULT_RULES_PATH = DEFAULT_RULES_PATH
    PRESET_DIR = PRESET_DIR
    USER_CONFIG_PATH = USER_RULES_PATH

    def __init__(self, config_path: Optional[Path] = None):
        self.config_path = config_path or self.USER_CONFIG_PATH
        self.rules = self.load_rules()

    def load_rules(self) -> RuleSet:
        """Load rules from user config or defaults (single source of truth)."""
        return load_active_ruleset(user_path=self.config_path, default_path=self.DEFAULT_RULES_PATH)

    def save_rules(self) -> None:
        """Save current rules to user config file."""
        save_ruleset(self.rules, self.config_path)

    def get_rules(self) -> RuleSet:
        return self.rules

    def update_rule(self, name: str, updates: Dict[str, Any]) -> None:
        """Update a specific rule's configuration."""
        if name in self.rules.rules:
            current = self.rules.rules[name]
            merged = current.model_dump()
            if "enabled" in updates:
                merged["enabled"] = updates["enabled"]
            if "weight" in updates:
                merged["weight"] = updates["weight"]
            if "parameters" in updates:
                params = merged.get("parameters", {})
                params.update(updates["parameters"])
                merged["parameters"] = params
            if "description" in updates:
                merged["description"] = updates["description"]
            if "ui" in updates:
                ui = merged.get("ui", {})
                ui.update(updates["ui"])
                merged["ui"] = ui
            self.rules.rules[name] = RuleParameter(**merged)
            self.save_rules()

    def list_presets(self) -> list[str]:
        return list_presets(self.PRESET_DIR)

    def apply_preset(self, name: str) -> RuleSet:
        preset_path = self.PRESET_DIR / f"{name}.json"
        preset_rules = load_ruleset(preset_path)
        if preset_rules is None:
            raise FileNotFoundError(f"Preset not found: {name}")
        self.rules = preset_rules
        self.save_rules()
        return self.rules

    def reset_to_default(self) -> RuleSet:
        default_rules = load_ruleset(self.DEFAULT_RULES_PATH)
        self.rules = default_rules or RuleSet()
        self.save_rules()
        return self.rules

    def import_rules(self, path: Path) -> RuleSet:
        rules = load_ruleset(path)
        if rules is None:
            raise ValueError(f"Invalid rules file: {path}")
        self.rules = rules
        self.save_rules()
        return self.rules
