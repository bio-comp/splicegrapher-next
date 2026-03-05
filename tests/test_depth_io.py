"""Parity tests for SGN depth record I/O helpers."""

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


def test_depth_io_is_depths_file_accepts_seekable_stream() -> None:
    from SpliceGrapher.formats.depth_io import is_depths_file

    stream = io.StringIO("C\tchr1\t3\nD\tchr1\t1:0,2:2\n")

    assert is_depths_file(stream) is True
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


def test_depth_io_read_depths_parses_junction_records(tmp_path: Path) -> None:
    from SpliceGrapher.formats.depth_io import read_depths
    from SpliceGrapher.formats.junction import SpliceJunction, parse_junction_record

    depths_path = tmp_path / "junction.depths"
    junction = SpliceJunction("chr1", 2, 5, (2, 3), "K", "+")
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

    depth_map, junction_map = read_depths(depths_path, parse_junction=parse_junction_record)

    assert set(depth_map) == {"chr1"}
    assert numpy.array_equal(
        depth_map["chr1"], numpy.array([0, 0, 2, 2, 2, 0, 0, 0], dtype=numpy.int32)
    )
    assert sorted(_junction_signature(j) for j in junction_map["chr1"]) == sorted(
        [_junction_signature(junction)]
    )


def test_depth_io_read_depths_respects_maxpos_without_extending_depth_array() -> None:
    from SpliceGrapher.formats.depth_io import read_depths
    from SpliceGrapher.formats.junction import parse_junction_record

    stream = io.StringIO("C\tchr1\t10\nD\tchr1\t10:5\n")

    depths, junctions = read_depths(
        stream,
        parse_junction=parse_junction_record,
        maxpos=3,
        junctions=False,
    )

    assert "chr1" in depths
    assert len(depths["chr1"]) == 3
    assert numpy.array_equal(depths["chr1"], numpy.array([5, 5, 5], dtype=numpy.int32))
    assert junctions == {}


def test_depth_io_requires_junction_parser_when_junctions_enabled() -> None:
    from SpliceGrapher.formats.depth_io import read_depths

    stream = io.StringIO("C\tchr1\t3\nD\tchr1\t1:0,2:2\n")

    try:
        read_depths(stream)
    except ValueError as exc:
        assert "parse_junction is required" in str(exc)
    else:
        raise AssertionError("Expected ValueError when parse_junction is missing")
