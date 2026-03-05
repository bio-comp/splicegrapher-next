"""Compatibility boundary for legacy ``shared.ShortRead`` APIs."""

from __future__ import annotations

import sys
from os import PathLike
from typing import BinaryIO, TextIO, TypeAlias

from SpliceGrapher.core.enums import JunctionCode, ShortReadCode
from SpliceGrapher.formats.depth_io import (
    DepthMap,
    DepthValues,
    is_depths_file,
)
from SpliceGrapher.formats.depth_io import (
    read_depths as _read_depths,
)
from SpliceGrapher.shared.ShortRead import SpliceJunction

JunctionMap: TypeAlias = dict[str, list[SpliceJunction]]
DepthSource: TypeAlias = str | PathLike[str] | TextIO | BinaryIO
_JUNCTION_CODES = frozenset(code.value for code in JunctionCode)


def _coerce_int(value: object, *, default: int) -> int:
    if value is None:
        return default
    if isinstance(value, int):
        return value
    if isinstance(value, str):
        return int(value)
    raise ValueError(f"Expected int-compatible value, got {type(value)!r}")


def _coerce_bool(value: object, *, default: bool) -> bool:
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, int):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "t", "yes", "y"}
    raise ValueError(f"Expected bool-compatible value, got {type(value)!r}")


def _parse_legacy_junction(record: str) -> SpliceJunction:
    """Parse one legacy ``J`` depth record into a ``SpliceJunction``."""
    parts = record.split("\t")
    if len(parts) != 9:
        raise ValueError(f"Invalid SpliceJunction record has {len(parts)} columns (expected 9)")

    if parts[0] != ShortReadCode.JUNCTION.value:
        raise ValueError(f"Invalid SpliceJunction code: {parts[0]}")

    if parts[7] not in _JUNCTION_CODES:
        raise ValueError(f"Invalid SpliceJunction type: {parts[7]}")

    chromosome = parts[1]
    strand = parts[2]
    positions = [int(value) for value in parts[3:7]]
    p1 = min(positions[0], positions[1])
    p2 = max(positions[0], positions[1])
    anchors = (positions[2], positions[3])

    junction = SpliceJunction(chromosome, p1, p2, anchors, parts[7], strand)
    junction.count = int(parts[8])
    return junction


def read_depths(source: DepthSource, **args: object) -> tuple[DepthMap, JunctionMap]:
    """Read SGN depth records while preserving legacy kwargs compatibility."""
    return _read_depths(
        source,
        parse_junction=_parse_legacy_junction,
        maxpos=_coerce_int(args.get("maxpos"), default=sys.maxsize),
        minanchor=_coerce_int(args.get("minanchor"), default=0),
        minjct=_coerce_int(args.get("minjct"), default=1),
        depths=_coerce_bool(args.get("depths"), default=True),
        junctions=_coerce_bool(args.get("junctions"), default=True),
        verbose=_coerce_bool(args.get("verbose"), default=False),
    )


__all__ = [
    "DepthMap",
    "DepthSource",
    "DepthValues",
    "JunctionMap",
    "SpliceJunction",
    "is_depths_file",
    "read_depths",
]
