"""Typed depth-record parsing helpers extracted from legacy ShortRead."""

from __future__ import annotations

import os
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

    count: int
    minpos: int

    def minAnchor(self) -> int: ...  # noqa: N802


TJunction = TypeVar("TJunction", bound=JunctionRecord)
JunctionMap: TypeAlias = dict[str, list[TJunction]]
ParseJunction = Callable[[str], TJunction]


def _as_text(line: str | bytes) -> str:
    return line.decode("utf-8") if isinstance(line, bytes) else line


def _load_lines(source: DepthSource) -> Iterator[str]:
    if isinstance(source, (str, os.PathLike)):
        with ez_open(Path(source)) as stream:
            for line in stream:
                yield _as_text(line)
        return
    if hasattr(source, "__iter__"):
        for line in source:
            yield _as_text(line)
        return
    if hasattr(source, "readline"):
        while True:
            line = source.readline()
            if line == "" or line == b"":
                break
            yield _as_text(line)
        return
    raise ValueError(f"Unrecognized depth record source: {type(source)!r}")


def is_depths_file(
    source: DepthSource,
    *,
    depth_codes: tuple[ShortReadCode, ...] | tuple[ShortReadCode, ShortReadCode, ShortReadCode] = (
        ShortReadCode.CHROM,
        ShortReadCode.DEPTH,
        ShortReadCode.JUNCTION,
    ),
) -> bool:
    """Return ``True`` when source begins with a valid SGN depth record code."""
    first_line: str | bytes
    if isinstance(source, (str, os.PathLike)):
        path = os.fspath(source)
        if not os.path.isfile(path):
            return False
        with ez_open(path) as stream:
            first_line = stream.readline()
    elif hasattr(source, "read"):
        reset_position = None
        if hasattr(source, "tell") and hasattr(source, "seek"):
            try:
                reset_position = source.tell()
            except (OSError, ValueError):
                reset_position = None
        first_line = source.readline()
        if reset_position is not None:
            source.seek(reset_position)
    else:
        return False

    parts = _as_text(first_line).strip().split("\t")
    return bool(parts) and parts[0] in depth_codes


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

        record_type = parts[0]
        chrom = parts[1]
        if record_type == chrom_code:
            if not depths:
                continue
            chromosome_limit = min(maxpos, int(parts[2]))
            chromosome_limits[chrom] = chromosome_limit
            depth_map[chrom] = numpy.zeros(chromosome_limit, dtype=numpy.int32)
            continue

        if record_type == jct_code:
            if not junctions:
                continue
            assert parse_junction is not None
            junction = parse_junction(text_line.strip())
            if junction.count < minjct:
                continue
            if junction.minAnchor() < minanchor:
                continue
            chromosome_limit = chromosome_limits.get(chrom, maxpos)
            if junction.minpos > chromosome_limit:
                continue
            junction_map.setdefault(chrom, []).append(junction)
            continue

        if not depths:
            continue
        if chrom not in depth_map:
            raise ValueError(f"No chromosome information specified for {chrom} depths")

        run_length_string = parts[2]
        position = 0
        chrom_limit = len(depth_map[chrom])
        for run_length_pair in run_length_string.split(","):
            run_length_str, height_str = run_length_pair.split(":", 1)
            run_length = int(run_length_str)
            height = int(height_str)
            upper_bound = min(chrom_limit, position + run_length)
            if upper_bound > position:
                depth_map[chrom][position:upper_bound] = height
            position += run_length
            if position >= chrom_limit:
                break

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
