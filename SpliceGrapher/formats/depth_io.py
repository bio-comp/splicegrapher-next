"""Typed depth-record parsing helpers extracted from legacy ShortRead."""

from __future__ import annotations

import sys
from collections.abc import Callable, Iterator
from os import PathLike
from pathlib import Path
from typing import BinaryIO, Protocol, TextIO, TypeAlias, TypeVar

import numpy

from SpliceGrapher.core.enums import ShortReadCode
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.progress import ProgressIndicator

DepthValues: TypeAlias = numpy.ndarray
DepthMap: TypeAlias = dict[str, DepthValues]
DepthSource: TypeAlias = str | PathLike[str] | TextIO | BinaryIO


class JunctionRecord(Protocol):
    """Minimum protocol needed to filter parsed junction records."""

    @property
    def count(self) -> int: ...

    @property
    def minpos(self) -> int: ...

    def minAnchor(self) -> int: ...  # noqa: N802


TJunction = TypeVar("TJunction", bound=JunctionRecord)
JunctionMap: TypeAlias = dict[str, list[TJunction]]
ParseJunction: TypeAlias = Callable[[str], TJunction]


def _as_text(line: str | bytes) -> str:
    return line.decode("utf-8") if isinstance(line, bytes) else line


def _iter_text_lines(stream: TextIO | BinaryIO) -> Iterator[str]:
    for line in stream:
        yield _as_text(line)


def _load_lines(source: DepthSource) -> Iterator[str]:
    """Yield decoded lines from a path or an already-open stream."""
    if isinstance(source, (str, PathLike)):
        with ez_open(Path(source)) as stream:
            yield from _iter_text_lines(stream)
        return

    yield from _iter_text_lines(source)


def is_depths_file(
    source: DepthSource,
    *,
    depth_codes: tuple[ShortReadCode, ...] = (
        ShortReadCode.CHROM,
        ShortReadCode.DEPTH,
        ShortReadCode.JUNCTION,
    ),
) -> bool:
    """Return ``True`` when source begins with a valid SGN depth record code."""
    first_line: str | bytes | None = None
    if isinstance(source, (str, PathLike)):
        path = Path(source)
        if not path.is_file():
            return False
        with ez_open(path) as stream:
            first_line = stream.readline()
    else:
        try:
            reset_position = source.tell()
            first_line = source.readline()
            source.seek(reset_position)
        except (AttributeError, OSError, ValueError):
            # For unseekable streams, avoid consuming input during format probes.
            # Callers should route those sources directly into read_depths().
            return False

    if not first_line:
        return False

    parts = _as_text(first_line).strip().split("\t")
    return bool(parts) and parts[0] in depth_codes


def _apply_run_length_depths(chrom_data: DepthValues, run_length_string: str) -> None:
    """Decode one run-length depth record into a chromosome depth array."""
    position = 0
    chrom_limit = len(chrom_data)
    for run_length_pair in run_length_string.split(","):
        run_length_str, height_str = run_length_pair.split(":", 1)
        run_length = int(run_length_str)
        height = int(height_str)

        upper_bound = min(chrom_limit, position + run_length)
        if upper_bound > position:
            chrom_data[position:upper_bound] = height
        position += run_length
        if position >= chrom_limit:
            break


def read_depths(
    source: DepthSource,
    *,
    parse_junction: ParseJunction[TJunction] | None = None,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    depths: bool = True,
    junctions: bool = True,
    verbose: bool = False,
    chrom_code: ShortReadCode = ShortReadCode.CHROM,
    jct_code: ShortReadCode = ShortReadCode.JUNCTION,
) -> tuple[DepthMap, JunctionMap[TJunction]]:
    """Read SGN depth records into per-chromosome depth and junction maps."""
    if junctions and parse_junction is None:
        raise ValueError("parse_junction is required when junctions=True")

    depth_map: DepthMap = {}
    junction_map: JunctionMap[TJunction] = {}
    chromosome_limits: dict[str, int] = {}

    indicator = ProgressIndicator(1000000, verbose=verbose)
    for text_line in _load_lines(source):
        indicator.update()
        parts = text_line.strip().split("\t")
        if len(parts) < 3:
            raise ValueError(f"Bad depths record at line {indicator.ctr}:\n{text_line}")

        record_type, chrom = parts[0], parts[1]
        if record_type == chrom_code:
            if depths:
                chromosome_limit = min(maxpos, int(parts[2]))
                chromosome_limits[chrom] = chromosome_limit
                depth_map[chrom] = numpy.zeros(chromosome_limit, dtype=numpy.int32)
            continue

        if record_type == jct_code:
            if junctions:
                assert parse_junction is not None
                junction = parse_junction(text_line.strip())
                if junction.count >= minjct and junction.minAnchor() >= minanchor:
                    chromosome_limit = chromosome_limits.get(chrom, maxpos)
                    if junction.minpos <= chromosome_limit:
                        junction_map.setdefault(chrom, []).append(junction)
            continue

        if not depths:
            continue
        if chrom not in depth_map:
            raise ValueError(f"No chromosome information specified for {chrom} depths")

        _apply_run_length_depths(depth_map[chrom], parts[2])

    indicator.finish()
    return depth_map, junction_map


__all__ = [
    "DepthMap",
    "DepthSource",
    "DepthValues",
    "JunctionMap",
    "JunctionRecord",
    "is_depths_file",
    "read_depths",
]
