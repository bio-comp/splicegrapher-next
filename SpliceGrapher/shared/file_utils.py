"""Filesystem and file-format helpers extracted from shared.utils."""

from __future__ import annotations

import gzip
from pathlib import Path
from typing import TextIO

_GZIP_MAGIC = b"\x1f\x8b"


def _as_path(path: str | Path) -> Path:
    return path if isinstance(path, Path) else Path(path)


def _is_gzip_file(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == _GZIP_MAGIC


def ez_open(file_name: str | Path) -> TextIO:
    """Open plain-text or gzipped files in text mode."""
    path = _as_path(file_name)
    if not path.is_file():
        raise FileNotFoundError(f"File does not exist at {path}")

    if path.suffix.lower() == ".gz":
        return gzip.open(path, mode="rt")

    if _is_gzip_file(path):
        return gzip.open(path, mode="rt")

    return path.open("r")


def file_len(path: str | Path) -> int:
    """Return an exact line count for a text file using buffered byte reads."""
    file_path = _as_path(path)
    line_count = 0
    saw_bytes = False
    last_chunk_ended_with_newline = True

    with file_path.open("rb") as instream:
        while chunk := instream.read(1024 * 1024):
            saw_bytes = True
            line_count += chunk.count(b"\n")
            last_chunk_ended_with_newline = chunk.endswith(b"\n")

    if saw_bytes and not last_chunk_ended_with_newline:
        line_count += 1
    return line_count


def file_prefix(path: str | Path) -> str:
    """Return the filename stem for a file path."""
    return _as_path(path).stem


def find_file(name: str, search_path: str, delim: str = ":") -> str | None:
    """Find the first matching file in a delimiter-separated path string."""
    for raw_path in search_path.split(delim):
        if not raw_path:
            continue
        base = Path(raw_path)
        file_path = base / name
        if file_path.is_file():
            return str(file_path)
    return None


def make_graph_list_file(splice_graph_dir: str | Path) -> str:
    """Write a sorted list of graph `.gff` paths under `<root>/*/*.gff`."""
    root = _as_path(splice_graph_dir)
    validate_dir(root)
    graph_list = sorted(str(path) for path in root.glob("*/*.gff") if path.is_file())
    if not graph_list:
        raise ValueError(f"No splice graphs found in {root}")

    result = Path(f"{root}.lis")
    result.write_text("\n".join(graph_list) + "\n")
    return str(result)


def validate_dir(path: str | Path) -> None:
    """Validate a directory path."""
    validate_file(path)
    if not _as_path(path).is_dir():
        raise NotADirectoryError(f"'{path}' is not a directory.")


def validate_file(path: str | Path) -> None:
    """Validate a non-empty existing path."""
    if not path:
        raise ValueError("Provided path is empty.")

    if not _as_path(path).exists():
        raise FileNotFoundError(f"File '{path}' not found.")


__all__ = [
    "ez_open",
    "file_len",
    "file_prefix",
    "find_file",
    "make_graph_list_file",
    "validate_dir",
    "validate_file",
]
