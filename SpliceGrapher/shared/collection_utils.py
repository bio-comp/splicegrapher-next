"""Collection and searching helpers extracted from shared.utils."""

from __future__ import annotations

import bisect
from collections.abc import Callable
from typing import Any, TypeVar

T = TypeVar("T")
CollectionValue = str | list[object] | set[object] | tuple[object, ...]


def as_list(value: CollectionValue, delim: str = ",") -> list[object]:
    """Create a list from a string, list, set, or tuple."""
    if isinstance(value, str):
        return value.split(delim)

    if isinstance(value, (list, set, tuple)):
        return list(value)

    raise TypeError(f"Expected a string, list, set, or tuple; received {type(value).__name__}")


def as_set(value: CollectionValue, delim: str = ",") -> set[object]:
    """Create a set from a string, list, set, or tuple."""
    if isinstance(value, str):
        return set(value.split(delim))

    if isinstance(value, (list, set, tuple)):
        return set(value)

    raise TypeError(f"Expected a string, list, set, or tuple; received {type(value).__name__}")


def binary_search(
    sequence: list[T],
    target: Any,
    key: Callable[[T], Any] | None = None,
    lo: int = 0,
    hi: int | None = None,
) -> int:
    """Return the bisect insertion index for ``target`` in ``sequence``."""
    if not sequence:
        raise ValueError("Cannot perform binary search on an empty sequence")

    if hi is None:
        hi = len(sequence)

    return bisect.bisect_left(sequence, target, lo=lo, hi=hi, key=key)


__all__ = [
    "as_list",
    "as_set",
    "binary_search",
]
