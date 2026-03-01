from __future__ import annotations

import tomllib
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
MANIFEST_PATH = REPO_ROOT / "tests" / "fixtures" / "parity_harness.toml"
FIXTURE_README_PATH = REPO_ROOT / "tests" / "fixtures" / "README.md"

REQUIRED_TEST_AREAS = {
    "core_shared_import_smoke",
    "annotation_io",
    "alignment_io_parity",
    "splicegrapher_alignment_io",
}


def test_parity_harness_manifest_declares_required_test_areas() -> None:
    manifest = tomllib.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    declared_areas = set(manifest["required_tests"].keys())
    missing = REQUIRED_TEST_AREAS - declared_areas
    assert not missing, f"Missing required parity test areas: {sorted(missing)}"


def test_parity_harness_manifest_points_to_existing_test_files() -> None:
    manifest = tomllib.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
    for area, relative_path in manifest["required_tests"].items():
        test_path = REPO_ROOT / relative_path
        assert test_path.exists(), f"{area} references missing file: {relative_path}"


def test_fixture_conventions_documentation_exists() -> None:
    assert FIXTURE_README_PATH.exists(), "Missing tests/fixtures/README.md conventions doc"
