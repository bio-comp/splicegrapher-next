"""Formatting and parsing helpers extracted from shared.utils."""

from __future__ import annotations

from collections.abc import Iterable
from datetime import datetime


def comma_format(value: int | float | str) -> str:
    """Format integer-compatible values using commas."""
    return f"{int(value):,}"


def dict_string(value_dict: dict[object, object], delim: str = ",") -> str:
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
    """Return a timestamp unique to the current second."""
    return datetime.now().strftime(format_string)


def time_string(message: str, format_string: str = "%X", trailing_newline: bool = False) -> str:
    """Return a message prefixed with a user-readable timestamp."""
    ts = datetime.now().strftime(format_string)
    result = f"{ts} {message}"
    if trailing_newline:
        result += "\n"
    return result


def to_numeric(value: str) -> int | float | str:
    """Attempt int then float conversion; return the original string on failure."""
    try:
        return int(value)
    except ValueError:
        pass

    try:
        return float(value)
    except ValueError:
        pass

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
