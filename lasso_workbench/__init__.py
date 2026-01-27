"""
Lasso Workbench - GUI application for lasso peptide discovery and dataset curation.

This package provides tools for:
- Predicting core sequences from precursor peptides
- ESM-2 semantic precursor discovery from GBK files
- Analyzing antiSMASH and MiBIG reports
- Tracking lab validation of predictions
- Managing and curating lasso peptide datasets
"""

__version__ = "0.3.0"
__author__ = "Magnus Ohle"

from pathlib import Path

# Package paths
PACKAGE_DIR = Path(__file__).parent
PROJECT_ROOT = PACKAGE_DIR.parent
DATA_DIR = PROJECT_ROOT / "data"
CONFIG_DIR = PACKAGE_DIR / "config"
