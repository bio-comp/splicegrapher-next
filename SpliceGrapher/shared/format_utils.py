"""Formatting and parsing helpers extracted from shared.utils."""

from __future__ import annotations

import re
from collections.abc import Iterable, Mapping
from datetime import datetime, timezone
from typing import TypeVar

_INT_PATTERN = re.compile(r"^[+-]?\d+$")
_FLOAT_PATTERN = re.compile(r"^[+-]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][+-]?\d+)?$")
_KeyT = TypeVar("_KeyT")
_ValueT = TypeVar("_ValueT")


def comma_format(value: int | float | str) -> str:
    """Format numeric values using thousands separators without truncation."""
    if isinstance(value, str):
        number = to_numeric(value)
        if isinstance(number, str):
            raise ValueError(f"Value is not numeric: {value!r}")
        return f"{number:,}"
    return f"{value:,}"


def dict_string(value_dict: Mapping[_KeyT, _ValueT], delim: str = ",") -> str:
    """Return a simple string representation for a mapping."""
    return delim.join(f"{key} -> {value}" for key, value in value_dict.items())


def list_string(values: Iterable[object], delim: str = ",") -> str:
    """Return a simple string representation for an iterable."""
    return delim.join(str(value) for value in values)


def substring_after(text: str, tag: str) -> str | None:
    """Return the substring after ``tag`` if present, else ``None``."""
    _, sep, after = text.partition(tag)
    return after if sep else None


def substring_before(text: str, tag: str) -> str | None:
    """Return the substring before ``tag`` if present, else ``None``."""
    before, sep, _ = text.partition(tag)
    return before if sep else None


def substring_between(text: str, tag1: str, tag2: str) -> str | None:
    """Return substring between ``tag1`` and ``tag2`` if both tags are present."""
    _, sep1, remainder = text.partition(tag1)
    if not sep1:
        return None

    between, sep2, _ = remainder.partition(tag2)
    return between if sep2 else None


def timestamp(format_string: str = "%Y%m%d%H%M%S") -> str:
    """Return a UTC timestamp unique to the current second."""
    return datetime.now(timezone.utc).strftime(format_string)


def time_string(
    message: str,
    format_string: str = "%H:%M:%S",
    trailing_newline: bool = False,
) -> str:
    """Return a message prefixed with a UTC timestamp."""
    ts = datetime.now(timezone.utc).strftime(format_string)
    result = f"{ts} {message}"
    if trailing_newline:
        result += "\n"
    return result


def to_numeric(value: str) -> int | float | str:
    """Attempt numeric conversion with a regex fast path, else return original string."""
    stripped = value.strip()
    if not stripped:
        return value

    if _INT_PATTERN.fullmatch(stripped):
        try:
            return int(stripped)
        except ValueError:
            return value

    if _FLOAT_PATTERN.fullmatch(stripped):
        try:
            return float(stripped)
        except ValueError:
            return value

    return value


__all__ = [
    "comma_format",
    "dict_string",
    "list_string",
    "substring_after",
    "substring_before",
    "substring_between",
    "timestamp",
    "time_string",
    "to_numeric",
]
