"""Compatibility boundary for legacy ``shared.ShortRead`` APIs.

This module centralizes legacy read-depth/junction imports so downstream
alignment logic can migrate away from ``shared.ShortRead`` incrementally.
"""

from __future__ import annotations

from os import PathLike
from typing import BinaryIO, TextIO, TypeAlias

from SpliceGrapher.shared.ShortRead import (
    SpliceJunction,
)
from SpliceGrapher.shared.ShortRead import (
    isDepthsFile as _is_depths_file_legacy,
)
from SpliceGrapher.shared.ShortRead import (
    readDepths as _read_depths_legacy,
)

DepthValues: TypeAlias = list[int]
DepthMap: TypeAlias = dict[str, DepthValues]
JunctionMap: TypeAlias = dict[str, list[SpliceJunction]]
DepthSource: TypeAlias = str | PathLike[str] | TextIO | BinaryIO


def is_depths_file(source: DepthSource) -> bool:
    """Return ``True`` when a source contains SGN depth records."""
    return bool(_is_depths_file_legacy(source))


def read_depths(source: DepthSource, **args: object) -> tuple[DepthMap, JunctionMap]:
    """Read SGN depth records through the legacy parser boundary."""
    depths, junctions = _read_depths_legacy(source, **args)
    return depths, junctions


def isDepthsFile(source: DepthSource) -> bool:
    """Compatibility wrapper for legacy camelCase API."""
    return is_depths_file(source)


def readDepths(source: DepthSource, **args: object) -> tuple[DepthMap, JunctionMap]:
    """Compatibility wrapper for legacy camelCase API."""
    return read_depths(source, **args)


__all__ = [
    "DepthMap",
    "DepthSource",
    "DepthValues",
    "JunctionMap",
    "SpliceJunction",
    "is_depths_file",
    "isDepthsFile",
    "read_depths",
    "readDepths",
]
