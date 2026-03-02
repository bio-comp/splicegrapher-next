"""Core shared types for SpliceGrapher."""

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AttrKey, EdgeType, NodeDisposition, RecordType, Strand
from SpliceGrapher.core.types import GenomicInterval, IntervalCoordinates

__all__ = [
    "AttrKey",
    "EdgeType",
    "GenomicInterval",
    "IntervalCoordinates",
    "NodeDisposition",
    "RecordType",
    "Strand",
    "coerce_enum",
]
