"""Smoke checks for core/shared extraction baseline."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def test_import_splicegraph_module() -> None:
    module = importlib.import_module("SpliceGrapher.SpliceGraph")
    assert hasattr(module, "SpliceGraph")


def test_import_shared_utils_and_config() -> None:
    utils_module = importlib.import_module("SpliceGrapher.shared.utils")
    config_module = importlib.import_module("SpliceGrapher.shared.config")

    assert hasattr(utils_module, "ProgressIndicator")
    assert hasattr(config_module, "load_config")
