from __future__ import annotations

import gzip
import importlib
from pathlib import Path

import pytest

from SpliceGrapher.shared.file_utils import (
    ez_open,
    file_len,
    file_prefix,
    find_file,
    make_graph_list_file,
    validate_dir,
    validate_file,
)


def test_ez_open_reads_plain_text(tmp_path: Path) -> None:
    plain = tmp_path / "example.txt"
    plain.write_text("a\nb\n")

    with ez_open(plain) as handle:
        assert handle.read() == "a\nb\n"


def test_ez_open_reads_gzip_text(tmp_path: Path) -> None:
    gz_path = tmp_path / "example.txt.gz"
    with gzip.open(gz_path, mode="wt") as handle:
        handle.write("x\ny\n")

    with ez_open(gz_path) as handle:
        assert handle.read() == "x\ny\n"


def test_ez_open_missing_file_raises(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        ez_open(tmp_path / "missing.txt")


def test_file_len_and_prefix(tmp_path: Path) -> None:
    file_path = tmp_path / "sample.data.txt"
    file_path.write_text("1\n2\n3\n")

    assert file_len(file_path) == 3
    assert file_prefix(file_path) == "sample.data"


def test_find_file(tmp_path: Path) -> None:
    one = tmp_path / "one"
    two = tmp_path / "two"
    one.mkdir()
    two.mkdir()

    target = two / "target.txt"
    target.write_text("ok")

    found = find_file("target.txt", f"{one}:{two}")
    assert found == str(target)

    assert find_file("missing.txt", f"{one}:{two}") is None


def test_make_graph_list_file(tmp_path: Path) -> None:
    root = tmp_path / "graphs"
    (root / "chr1").mkdir(parents=True)
    (root / "chr2").mkdir(parents=True)

    g1 = root / "chr1" / "a.gff"
    g2 = root / "chr2" / "b.gff"
    g1.write_text("g1\n")
    g2.write_text("g2\n")

    list_path = make_graph_list_file(root)
    result = Path(list_path)
    assert result.exists()

    lines = result.read_text().splitlines()
    assert lines == sorted([str(g1), str(g2)])


def test_validate_file_and_dir(tmp_path: Path) -> None:
    data_file = tmp_path / "file.txt"
    data_file.write_text("data")

    validate_file(data_file)
    validate_dir(tmp_path)

    with pytest.raises(ValueError):
        validate_file("")

    with pytest.raises(FileNotFoundError):
        validate_file(tmp_path / "does_not_exist")

    with pytest.raises(NotADirectoryError):
        validate_dir(data_file)


def test_legacy_wrappers_are_not_exposed() -> None:
    module = importlib.import_module("SpliceGrapher.shared.file_utils")
    legacy_names = [
        "ezopen",
        "fileLen",
        "filePrefix",
        "findFile",
        "makeGraphListFile",
        "validateDir",
        "validateFile",
    ]
    for name in legacy_names:
        assert not hasattr(module, name)
