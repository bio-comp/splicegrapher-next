"""Parity tests for extracted SGN depth record I/O helpers."""

from __future__ import annotations

import io
from pathlib import Path


def _junction_signature(junction: object) -> tuple[object, ...]:
    return (
        getattr(junction, "chromosome"),
        getattr(junction, "minpos"),
        getattr(junction, "maxpos"),
        tuple(getattr(junction, "anchors")),
        getattr(junction, "sjCode"),
        getattr(junction, "count"),
        getattr(junction, "strand"),
    )


def test_depth_io_is_depths_file_matches_shortread_behavior() -> None:
    from SpliceGrapher.formats.depth_io import is_depths_file
    from SpliceGrapher.shared.ShortRead import isDepthsFile

    stream = io.StringIO("C\tchr1\t3\nD\tchr1\t1:0,2:2\n")

    assert is_depths_file(stream) is True
    assert isDepthsFile(stream) is True
    assert stream.tell() == 0


def test_depth_io_read_depths_matches_legacy_parser(tmp_path: Path) -> None:
    from SpliceGrapher.formats.depth_io import read_depths
    from SpliceGrapher.shared.ShortRead import (
        KNOWN_JCT,
        SpliceJunction,
        readDepths,
        stringToJunction,
    )

    depths_path = tmp_path / "junction.depths"
    junction = SpliceJunction("chr1", 2, 5, (2, 3), KNOWN_JCT, "+")
    junction.count = 7
    depths_path.write_text(
        "\n".join(
            [
                "C\tchr1\t8",
                "D\tchr1\t2:0,3:2,3:0",
                junction.toString(),
                "",
            ]
        ),
        encoding="utf-8",
    )

    legacy_depths, legacy_junctions = readDepths(depths_path)
    new_depths, new_junctions = read_depths(depths_path, parse_junction=stringToJunction)

    assert new_depths == legacy_depths
    assert sorted(_junction_signature(j) for j in new_junctions["chr1"]) == sorted(
        _junction_signature(j) for j in legacy_junctions["chr1"]
    )
