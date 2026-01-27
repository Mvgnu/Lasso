"""
Rule configuration I/O helpers.
Single source of truth: JSON on disk (defaults or user custom).
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

from lasso_workbench.schemas.config import RuleSet

DEFAULT_RULES_PATH = Path(__file__).resolve().parents[1] / "config" / "defaults" / "rules.json"
USER_RULES_PATH = Path.home() / ".lasso_workbench" / "custom_rules.json"
PRESET_DIR = Path(__file__).resolve().parents[1] / "config" / "presets"


def load_ruleset(path: Path) -> Optional[RuleSet]:
    if not path.exists():
        return None
    try:
        data = json.loads(path.read_text())
        return RuleSet(**data)
    except Exception:
        return None


def load_active_ruleset(
    user_path: Path = USER_RULES_PATH,
    default_path: Path = DEFAULT_RULES_PATH,
) -> RuleSet:
    user_rules = load_ruleset(user_path)
    if user_rules is not None:
        return user_rules
    default_rules = load_ruleset(default_path)
    return default_rules or RuleSet()


def save_ruleset(rules: RuleSet, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(rules.model_dump_json(indent=2))


def list_presets(preset_dir: Path = PRESET_DIR) -> list[str]:
    if not preset_dir.exists():
        return []
    return sorted(p.stem for p in preset_dir.glob("*.json"))
