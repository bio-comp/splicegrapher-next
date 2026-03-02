"""Compatibility boundary for legacy ``shared.ShortRead`` APIs."""

from __future__ import annotations

import sys
from os import PathLike
from typing import BinaryIO, TextIO, TypeAlias

from SpliceGrapher.formats.depth_io import (
    DepthMap,
    DepthValues,
    is_depths_file,
)
from SpliceGrapher.formats.depth_io import (
    read_depths as _read_depths,
)
from SpliceGrapher.shared.ShortRead import SpliceJunction, stringToJunction

JunctionMap: TypeAlias = dict[str, list[SpliceJunction]]
DepthSource: TypeAlias = str | PathLike[str] | TextIO | BinaryIO


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


def read_depths(source: DepthSource, **args: object) -> tuple[DepthMap, JunctionMap]:
    """Read SGN depth records while preserving legacy kwargs compatibility."""
    return _read_depths(
        source,
        parse_junction=stringToJunction,
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
