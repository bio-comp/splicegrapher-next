"""Shared type aliases and protocols for alignment I/O."""

from __future__ import annotations

import io
from os import PathLike
from typing import Protocol, TypeAlias

import numpy

from SpliceGrapher.formats.depth_io import DepthSource
from SpliceGrapher.formats.junction import SpliceJunction

ChromosomeInput: TypeAlias = str | list[str] | tuple[str, ...] | set[str] | None
ReferencePath: TypeAlias = str | PathLike[str] | None
SamLine: TypeAlias = str | bytes
SamLineCollection: TypeAlias = list[SamLine] | tuple[SamLine, ...] | set[SamLine]
SamTextSource: TypeAlias = io.IOBase | SamLineCollection
AlignmentPath: TypeAlias = str | PathLike[str]
AlignmentSource: TypeAlias = AlignmentPath | SamTextSource
ReadDataSource: TypeAlias = AlignmentSource | DepthSource
DepthMap: TypeAlias = dict[str, numpy.ndarray]
JunctionMap: TypeAlias = dict[str, list[SpliceJunction]]
AlignmentMap: TypeAlias = dict[str, list[tuple[int, int]]]
JunctionKey: TypeAlias = tuple[str, str, int, int]
CollectResult: TypeAlias = tuple[DepthMap, JunctionMap]
CollectResultWithAlignments: TypeAlias = tuple[DepthMap, JunctionMap, AlignmentMap]


class GeneBounds(Protocol):
    id: str
    strand: str
    minpos: int
    maxpos: int
