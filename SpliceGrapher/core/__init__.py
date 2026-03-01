"""Core shared types for SpliceGrapher."""

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AttrKey, EdgeType, NodeDisposition, RecordType, Strand

__all__ = [
    "AttrKey",
    "EdgeType",
    "NodeDisposition",
    "RecordType",
    "Strand",
    "coerce_enum",
]
