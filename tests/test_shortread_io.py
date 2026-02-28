"""Parity-focused tests for shared.ShortRead modernization."""

from __future__ import annotations

import io
from pathlib import Path

from SpliceGrapher.shared.ShortRead import isDepthsFile, writeDepths


def test_shortread_source_does_not_use_fasta_star_import() -> None:
    shortread_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "shared" / "ShortRead.py"
    )
    source = shortread_path.read_text(encoding="utf-8")

    assert "from SpliceGrapher.formats.fasta import *" not in source


def test_is_depths_file_accepts_file_like_stream() -> None:
    stream = io.StringIO("C\tchr1\t3\nD\tchr1\t1:0,2:2\n")

    assert isDepthsFile(stream) is True
    assert stream.tell() == 0


def test_write_depths_accepts_text_io_stream() -> None:
    out_stream = io.StringIO()

    writeDepths(out_stream, {"chr1": [0, 2, 2]}, verbose=False)
    payload = out_stream.getvalue()

    assert payload.startswith("C\tchr1\t3\n")
    assert "D\tchr1\t" in payload
