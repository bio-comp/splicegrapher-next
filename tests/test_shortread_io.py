"""Parity-focused tests for shared.ShortRead modernization."""

from __future__ import annotations

import inspect
import io
from collections.abc import Callable
from pathlib import Path
from typing import cast

import pytest

from SpliceGrapher.shared import ShortRead as shortread
from SpliceGrapher.shared.ShortRead import Cluster, Read, SpliceJunction, isDepthsFile, writeDepths


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


def test_shortread_source_no_longer_uses_python2_object_bases() -> None:
    shortread_path = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "shared" / "ShortRead.py"
    )
    source = shortread_path.read_text(encoding="utf-8")

    assert "class Read(object):" not in source
    assert "class ReadPair(object):" not in source
    assert "class Cluster(object):" not in source


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


def test_write_depths_uses_non_mutable_default_for_junction_map() -> None:
    signature = inspect.signature(shortread.writeDepths)
    assert signature.parameters["jctDict"].default is None


def test_shortread_read_objects_sort_by_chromosome_and_coordinates() -> None:
    reads = [
        Read("chr2", 10, 12, "+"),
        Read("chr1", 20, 22, "+"),
        Read("chr1", 5, 7, "+"),
    ]

    ordered = sorted(reads)

    assert [r.chromosome for r in ordered] == ["chr1", "chr1", "chr2"]
    assert [r.minpos for r in ordered] == [5, 20, 10]


def test_shortread_cluster_objects_sort_by_coordinates() -> None:
    cluster_a = Cluster("chr1", 20, 1, id="a")
    cluster_a.addPosition(22, 1)
    cluster_b = Cluster("chr1", 5, 1, id="b")
    cluster_b.addPosition(6, 1)

    ordered = sorted([cluster_a, cluster_b])

    assert [c.id for c in ordered] == ["b", "a"]


def test_write_depths_sorts_junctions_without_python2_cmp() -> None:
    out_stream = io.StringIO()
    jct_late = SpliceJunction("chr1", 30, 40, (2, 2), shortread.KNOWN_JCT, "+")
    jct_early = SpliceJunction("chr1", 10, 20, (2, 2), shortread.KNOWN_JCT, "+")

    writeDepths(
        out_stream,
        {"chr1": [0, 0, 1, 1]},
        {"chr1": [jct_late, jct_early]},
        verbose=False,
    )

    assert [j.minpos for j in [jct_late, jct_early]] == [30, 10]
    junction_lines = [line for line in out_stream.getvalue().splitlines() if line.startswith("J\t")]
    assert len(junction_lines) == 2
    assert int(junction_lines[0].split("\t")[3]) < int(junction_lines[1].split("\t")[3])
