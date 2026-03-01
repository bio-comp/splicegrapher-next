"""Collection and searching helpers extracted from shared.utils."""

from __future__ import annotations

import bisect
import warnings
from collections.abc import Callable
from typing import Any, TypeVar

from SpliceGrapher.shared.process_utils import getAttribute

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
    """Return the index in ``sequence`` whose key is closest to ``target``."""
    if not sequence:
        raise ValueError("Cannot perform binary search on an empty sequence")

    if hi is None:
        hi = len(sequence)

    idx = bisect.bisect_left(sequence, target, lo=lo, hi=hi, key=key)
    return min(idx, len(sequence) - 1)


def _warn_deprecated(old_name: str, new_name: str) -> None:
    warnings.warn(
        f"'{old_name}' is deprecated and will be removed in a future release. "
        f"Use '{new_name}' instead.",
        DeprecationWarning,
        stacklevel=3,
    )


# Compatibility wrappers retained while older modules migrate to snake_case names.
def asList(value: CollectionValue, delim: str = ",") -> list[object]:
    _warn_deprecated("asList", "as_list")
    return as_list(value, delim)


def asSet(value: CollectionValue, delim: str = ",") -> set[object]:
    _warn_deprecated("asSet", "as_set")
    return as_set(value, delim)


def bsearch(
    sequence: list[T],
    target: Any,
    getValue: Callable[[T], Any] = lambda item: item,
    **args: Any,
) -> int:
    _warn_deprecated("bsearch", "binary_search")
    lo = getAttribute("low", 0, **args)
    high = getAttribute("high", len(sequence) - 1, **args)
    return binary_search(sequence, target, key=getValue, lo=lo, hi=high + 1)


__all__ = [
    "as_list",
    "as_set",
    "binary_search",
    "asList",
    "asSet",
    "bsearch",
]
