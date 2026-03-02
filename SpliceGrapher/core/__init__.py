"""Core shared types for SpliceGrapher."""

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AttrKey, EdgeType, NodeDisposition, RecordType, Strand
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    IntervalBounds,
    batch_overlaps,
    interval_contains,
    intervals_overlap,
)
from SpliceGrapher.core.types import GenomicInterval, IntervalCoordinates

__all__ = [
    "AttrKey",
    "EdgeType",
    "GenomicInterval",
    "InMemoryIntervalIndex",
    "IntervalBounds",
    "IntervalCoordinates",
    "NodeDisposition",
    "RecordType",
    "Strand",
    "batch_overlaps",
    "coerce_enum",
    "interval_contains",
    "intervals_overlap",
]
