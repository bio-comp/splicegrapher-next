"""Checks for direct focused shared-module usage."""

from __future__ import annotations

import importlib


def test_shared_modules_export_expected_helpers() -> None:
    file_utils_module = importlib.import_module("SpliceGrapher.shared.file_utils")
    format_utils_module = importlib.import_module("SpliceGrapher.shared.format_utils")
    collection_utils_module = importlib.import_module("SpliceGrapher.shared.collection_utils")
    process_module = importlib.import_module("SpliceGrapher.shared.process")
    progress_module = importlib.import_module("SpliceGrapher.shared.progress")

    assert hasattr(file_utils_module, "ez_open")
    assert hasattr(format_utils_module, "comma_format")
    assert hasattr(collection_utils_module, "as_list")
    assert hasattr(process_module, "get_attribute")
    assert hasattr(progress_module, "ProgressIndicator")


def test_fasta_readers_use_focused_file_utils_module() -> None:
    fasta_module = importlib.import_module("SpliceGrapher.formats.fasta")
    readers_module = importlib.import_module("SpliceGrapher.formats.fasta.readers")
    file_utils_module = importlib.import_module("SpliceGrapher.shared.file_utils")

    assert not hasattr(fasta_module, "ez_open")
    assert readers_module.ez_open is file_utils_module.ez_open


def test_representative_focused_utility_behavior() -> None:
    collection_utils_module = importlib.import_module("SpliceGrapher.shared.collection_utils")
    format_utils_module = importlib.import_module("SpliceGrapher.shared.format_utils")

    assert collection_utils_module.as_list("a,b,c") == ["a", "b", "c"]
    assert collection_utils_module.as_set("a,b,a") == {"a", "b"}
    assert format_utils_module.to_numeric("12") == 12
    assert format_utils_module.to_numeric("12.5") == 12.5
