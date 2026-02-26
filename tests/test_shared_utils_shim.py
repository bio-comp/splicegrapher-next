"""Compatibility checks for the shared.utils decomposition shim."""

from __future__ import annotations

import importlib


def test_shim_reexports_match_focused_modules() -> None:
    utils_module = importlib.import_module("SpliceGrapher.shared.utils")
    file_utils_module = importlib.import_module("SpliceGrapher.shared.file_utils")
    format_utils_module = importlib.import_module("SpliceGrapher.shared.format_utils")
    process_utils_module = importlib.import_module("SpliceGrapher.shared.process_utils")
    progress_module = importlib.import_module("SpliceGrapher.shared.progress")

    assert utils_module.ezopen is file_utils_module.ezopen
    assert utils_module.commaFormat is format_utils_module.commaFormat
    assert utils_module.getAttribute is process_utils_module.getAttribute
    assert utils_module.ProgressIndicator is progress_module.ProgressIndicator


def test_fasta_import_uses_focused_file_utils_module() -> None:
    fasta_module = importlib.import_module("SpliceGrapher.formats.fasta")
    file_utils_module = importlib.import_module("SpliceGrapher.shared.file_utils")

    assert fasta_module.ezopen is file_utils_module.ezopen


def test_shim_retains_representative_utility_behavior() -> None:
    utils_module = importlib.import_module("SpliceGrapher.shared.utils")

    assert utils_module.asList("a,b,c") == ["a", "b", "c"]
    assert utils_module.asSet("a,b,a") == {"a", "b"}
    assert utils_module.to_numeric("12") == 12
    assert utils_module.to_numeric("12.5") == 12.5
