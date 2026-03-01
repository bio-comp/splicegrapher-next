"""Helpers for coercing external values into canonical enum domains."""

from __future__ import annotations

from enum import Enum
from typing import TypeVar

from SpliceGrapher.core.enums import AttrKey, EdgeType, NodeDisposition, RecordType, Strand

_EnumT = TypeVar("_EnumT", bound=Enum)


def _coerce_enum(value: str | _EnumT, enum_type: type[_EnumT], *, field: str) -> _EnumT:
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(value)
    except ValueError as exc:
        valid_values = ", ".join(item.value for item in enum_type)
        raise ValueError(
            f"Invalid {field} value {value!r}. Expected one of: {valid_values}"
        ) from exc


def coerce_strand(value: str | Strand) -> Strand:
    return _coerce_enum(value, Strand, field="strand")


def coerce_record_type(value: str | RecordType) -> RecordType:
    return _coerce_enum(value, RecordType, field="record_type")


def coerce_node_disposition(value: str | NodeDisposition) -> NodeDisposition:
    return _coerce_enum(value, NodeDisposition, field="node_disposition")


def coerce_edge_type(value: str | EdgeType) -> EdgeType:
    return _coerce_enum(value, EdgeType, field="edge_type")


def coerce_attr_key(value: str | AttrKey) -> AttrKey:
    return _coerce_enum(value, AttrKey, field="attribute_key")
