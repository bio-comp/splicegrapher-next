"""Helpers for coercing external values into canonical enum domains."""

from __future__ import annotations

from enum import Enum
from typing import TypeVar

_EnumT = TypeVar("_EnumT", bound=Enum)


def coerce_enum(value: str | _EnumT, enum_type: type[_EnumT], *, field: str) -> _EnumT:
    """Coerce an external value into a canonical enum domain."""
    if isinstance(value, enum_type):
        return value
    try:
        return enum_type(value)
    except ValueError as exc:
        valid_values = ", ".join(item.value for item in enum_type)
        raise ValueError(
            f"Invalid {field} value {value!r}. Expected one of: {valid_values}"
        ) from exc
