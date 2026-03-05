"""Parity-focused tests for shared.ShortRead modernization."""

from __future__ import annotations

import inspect
import io
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest

from SpliceGrapher.shared import ShortRead as shortread
from SpliceGrapher.shared.ShortRead import isDepthsFile, writeDepths


def test_shortread_source_does_not_use_fasta_star_import() -> None:
    shortread_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "shared" / "ShortRead.py"
    )
    source = shortread_path.read_text(encoding="utf-8")

    assert "from SpliceGrapher.formats.fasta import *" not in source


def test_shortread_source_no_longer_uses_process_utils_getattribute() -> None:
    shortread_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "shared" / "ShortRead.py"
    )
    source = shortread_path.read_text(encoding="utf-8")

    assert "from SpliceGrapher.shared.process_utils import getAttribute" not in source
    assert "getAttribute(" not in source


@pytest.mark.parametrize("func_name", ["depthsToClusters", "readDepths"])
def test_shortread_hot_path_signatures_are_explicit(func_name: str) -> None:
    signature = inspect.signature(getattr(shortread, func_name))
    assert all(
        param.kind != inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()
    )


def test_read_depths_rejects_unknown_keyword_argument() -> None:
    read_depths = cast(Callable[..., object], shortread.readDepths)
    with pytest.raises(TypeError):
        read_depths(io.StringIO("C\tchr1\t3\n"), nonsense=True)


def test_depths_to_clusters_rejects_unknown_keyword_argument() -> None:
    depths_to_clusters = cast(Callable[..., object], shortread.depthsToClusters)
    with pytest.raises(TypeError):
        depths_to_clusters("chr1", [0, 2, 2], nonsense=True)


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
