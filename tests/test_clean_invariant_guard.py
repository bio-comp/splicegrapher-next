"""Tests for SGN clean-invariant ratchet guard logic."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "ci" / "check_clean_invariant.py"


def _load_guard_module() -> ModuleType:
    if not SCRIPT_PATH.exists():
        raise RuntimeError("Expected scripts/ci/check_clean_invariant.py to exist")

    spec = importlib.util.spec_from_file_location("check_clean_invariant", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load clean invariant guard module")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_detect_ignore_growth_rejects_new_file_entries() -> None:
    module = _load_guard_module()

    baseline = {"a.py": {"E"}}
    current = {"a.py": {"E"}, "new.py": {"F"}}

    violations = module.detect_ignore_growth(current=current, baseline=baseline)

    assert any("new.py" in entry for entry in violations)


def test_detect_ignore_growth_rejects_widened_rule_sets() -> None:
    module = _load_guard_module()

    baseline = {"a.py": {"E"}}
    current = {"a.py": {"E", "F"}}

    violations = module.detect_ignore_growth(current=current, baseline=baseline)

    assert any("a.py" in entry and "F" in entry for entry in violations)


def test_detect_ignore_growth_allows_debt_reduction() -> None:
    module = _load_guard_module()

    baseline = {"a.py": {"E", "F"}}
    current = {"a.py": {"E"}}

    assert module.detect_ignore_growth(current=current, baseline=baseline) == []


def test_find_unexpected_executable_python_paths_flags_library_modules() -> None:
    module = _load_guard_module()

    ls_lines = [
        "100644 abc 0\tSpliceGrapher/shared/ok.py",
        "100755 abc 0\tSpliceGrapher/formats/fasta.py",
        "100755 abc 0\tscripts/release.py",
    ]

    unexpected = module.find_unexpected_executable_python_paths(
        ls_files_lines=ls_lines,
        executable_allowlist={"scripts/release.py"},
    )

    assert unexpected == ["SpliceGrapher/formats/fasta.py"]
