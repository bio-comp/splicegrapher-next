"""Retirement-boundary tests for legacy ShortRead compatibility modules."""

from __future__ import annotations

from pathlib import Path


def test_alignment_io_imports_modern_depth_and_junction_boundaries() -> None:
    package_dir = Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "alignment_io"
    source = "\n".join(
        path.read_text(encoding="utf-8") for path in sorted(package_dir.glob("*.py"))
    )

    assert "from SpliceGrapher.shared.ShortRead import" not in source
    assert "from SpliceGrapher.formats.shortread_compat import" not in source
    assert "from SpliceGrapher.formats.depth_io import" in source
    assert "from SpliceGrapher.formats.junction import" in source


def test_legacy_shortread_and_compat_modules_are_removed() -> None:
    repo_root = Path(__file__).resolve().parents[1]
    assert not (repo_root / "SpliceGrapher" / "formats" / "shortread_compat.py").exists()
    assert not (repo_root / "SpliceGrapher" / "shared" / "ShortRead.py").exists()
