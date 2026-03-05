"""Compatibility-boundary tests for ShortRead deprecation slices."""

from __future__ import annotations

from pathlib import Path


def test_alignment_io_imports_shortread_compat_boundary() -> None:
    alignment_io_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "alignment_io.py"
    )
    source = alignment_io_path.read_text(encoding="utf-8")

    assert "from SpliceGrapher.shared.ShortRead import" not in source
    assert "from SpliceGrapher.formats.shortread_compat import" in source


def test_shortread_compat_round_trips_depth_file(tmp_path: Path) -> None:
    from SpliceGrapher.formats.shortread_compat import is_depths_file, read_depths

    depths_path = tmp_path / "sample.depths"
    depths_path.write_text(
        "\n".join(
            [
                "C\tchr1\t5",
                "D\tchr1\t2:0,3:4",
                "",
            ]
        ),
        encoding="utf-8",
    )

    assert is_depths_file(depths_path)
    depths, junctions = read_depths(depths_path)
    assert "chr1" in depths
    assert len(depths["chr1"]) == 5
    assert junctions == {}


def test_shortread_compat_parses_legacy_junction_records(tmp_path: Path) -> None:
    from SpliceGrapher.formats.shortread_compat import read_depths

    depths_path = tmp_path / "with_junction.depths"
    depths_path.write_text(
        "\n".join(
            [
                "C\tchr1\t8",
                "D\tchr1\t2:0,6:1",
                "J\tchr1\t+\t2\t5\t2\t3\tK\t7",
                "",
            ]
        ),
        encoding="utf-8",
    )

    _, junctions = read_depths(depths_path)
    assert "chr1" in junctions
    assert len(junctions["chr1"]) == 1
    junction = junctions["chr1"][0]
    assert junction.chromosome == "chr1"
    assert junction.minpos == 2
    assert junction.maxpos == 5
    assert tuple(junction.anchors) == (2, 3)
    assert junction.sjCode == "K"
    assert junction.count == 7
    assert junction.strand == "+"


def test_shortread_compat_imports_depth_io_boundary() -> None:
    compat_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "shortread_compat.py"
    )
    source = compat_path.read_text(encoding="utf-8")

    assert "from SpliceGrapher.formats.depth_io import" in source
    assert "stringToJunction" not in source
    assert "def isDepthsFile(" not in source
    assert "def readDepths(" not in source
