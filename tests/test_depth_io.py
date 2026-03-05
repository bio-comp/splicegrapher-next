"""Parity tests for extracted SGN depth record I/O helpers."""

from __future__ import annotations

import io
from pathlib import Path
from typing import TextIO, cast

import numpy


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


def test_depth_io_is_depths_file_rejects_unseekable_stream_without_consuming() -> None:
    from SpliceGrapher.formats.depth_io import is_depths_file

    class UnseekableStream:
        def __init__(self, text: str) -> None:
            self._stream = io.StringIO(text)

        def readline(self, size: int = -1) -> str:
            return self._stream.readline(size)

        def tell(self) -> int:
            msg = "stream is unseekable"
            raise OSError(msg)

        def seek(self, offset: int, whence: int = 0) -> int:
            msg = "stream is unseekable"
            raise OSError(msg)

    stream = UnseekableStream("C\tchr1\t3\nD\tchr1\t1:0,2:2\n")

    assert is_depths_file(cast(TextIO, stream)) is False
    assert stream.readline() == "C\tchr1\t3\n"


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

    assert new_depths.keys() == legacy_depths.keys()
    for chrom in new_depths:
        assert numpy.array_equal(new_depths[chrom], numpy.asarray(legacy_depths[chrom]))
    assert sorted(_junction_signature(j) for j in new_junctions["chr1"]) == sorted(
        _junction_signature(j) for j in legacy_junctions["chr1"]
    )


def test_depth_io_read_depths_respects_maxpos_without_extending_depth_array() -> None:
    from SpliceGrapher.formats.depth_io import read_depths
    from SpliceGrapher.shared.ShortRead import stringToJunction

    stream = io.StringIO("C\tchr1\t10\nD\tchr1\t10:5\n")

    depths, junctions = read_depths(
        stream,
        maxpos=3,
        junctions=False,
        parse_junction=stringToJunction,
    )

    assert "chr1" in depths
    assert len(depths["chr1"]) == 3
    assert numpy.array_equal(depths["chr1"], numpy.array([5, 5, 5], dtype=numpy.int32))
    assert junctions == {}
