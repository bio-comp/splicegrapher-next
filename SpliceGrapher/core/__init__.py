"""Core shared types for SpliceGrapher."""

from SpliceGrapher.core.enum_coercion import (
    coerce_attr_key,
    coerce_edge_type,
    coerce_node_disposition,
    coerce_record_type,
    coerce_strand,
)
from SpliceGrapher.core.enums import AttrKey, EdgeType, NodeDisposition, RecordType, Strand

__all__ = [
    "AttrKey",
    "EdgeType",
    "NodeDisposition",
    "RecordType",
    "Strand",
    "coerce_attr_key",
    "coerce_edge_type",
    "coerce_node_disposition",
    "coerce_record_type",
    "coerce_strand",
]
