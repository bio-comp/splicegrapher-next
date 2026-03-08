"""Smoke checks for core/shared extraction baseline."""

from __future__ import annotations

import importlib
import sys
from pathlib import Path

import pytest

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))


def test_import_splice_graph_core_module() -> None:
    module = importlib.import_module("SpliceGrapher.core.splice_graph")
    assert hasattr(module, "SpliceGraph")


def test_import_splice_graph_parser_module() -> None:
    module = importlib.import_module("SpliceGrapher.formats.parsers.splice_graph")
    assert hasattr(module, "SpliceGraphParser")


def test_import_shared_modules_without_utils_shim() -> None:
    file_utils_module = importlib.import_module("SpliceGrapher.shared.file_utils")
    process_utils_module = importlib.import_module("SpliceGrapher.shared.process_utils")
    progress_module = importlib.import_module("SpliceGrapher.shared.progress")
    config_module = importlib.import_module("SpliceGrapher.shared.config")

    assert hasattr(file_utils_module, "ez_open")
    assert hasattr(process_utils_module, "runCommand")
    assert hasattr(progress_module, "ProgressIndicator")
    assert hasattr(config_module, "load_config")


def test_shared_utils_shim_removed() -> None:
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("SpliceGrapher.shared.utils")
