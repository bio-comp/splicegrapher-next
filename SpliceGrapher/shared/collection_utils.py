"""Collection and searching helpers extracted from shared.utils."""

from __future__ import annotations

import bisect
from collections.abc import Callable
from typing import Protocol, TypeVar, cast, overload


class SupportsRichComparison(Protocol):
    def __lt__(self, other: object, /) -> bool: ...
    def __gt__(self, other: object, /) -> bool: ...


T = TypeVar("T")
K = TypeVar("K")
CollectionValue = str | list[T] | set[T] | tuple[T, ...]


@overload
def as_list(value: str, delim: str = ",") -> list[str]: ...


@overload
def as_list(value: list[T] | set[T] | tuple[T, ...], delim: str = ",") -> list[T]: ...


def as_list(value: CollectionValue[T], delim: str = ",") -> list[str] | list[T]:
    """Create a list from a string, list, set, or tuple."""
    if isinstance(value, str):
        return value.split(delim)

    if isinstance(value, (list, set, tuple)):
        return list(value)

    raise TypeError(f"Expected a string, list, set, or tuple; received {type(value).__name__}")


@overload
def as_set(value: str, delim: str = ",") -> set[str]: ...


@overload
def as_set(value: list[T] | set[T] | tuple[T, ...], delim: str = ",") -> set[T]: ...


def as_set(value: CollectionValue[T], delim: str = ",") -> set[str] | set[T]:
    """Create a set from a string, list, set, or tuple."""
    if isinstance(value, str):
        return set(value.split(delim))

    if isinstance(value, (list, set, tuple)):
        return set(value)

    raise TypeError(f"Expected a string, list, set, or tuple; received {type(value).__name__}")


def binary_search(
    sequence: list[T],
    target: K,
    key: Callable[[T], K] | None = None,
    lo: int = 0,
    hi: int | None = None,
) -> int:
    """Return the bisect insertion index for ``target`` in ``sequence``."""
    if not sequence:
        raise ValueError("Cannot perform binary search on an empty sequence")

    if hi is None:
        hi = len(sequence)

    comparable_target = cast(SupportsRichComparison, target)

    if key is None:
        ordered_sequence = cast(list[SupportsRichComparison], sequence)
        return bisect.bisect_left(ordered_sequence, comparable_target, lo=lo, hi=hi)

    keyed_sequence = [key(item) for item in sequence]
    comparable_keys = cast(list[SupportsRichComparison], keyed_sequence)
    return bisect.bisect_left(comparable_keys, comparable_target, lo=lo, hi=hi)


__all__ = [
    "as_list",
    "as_set",
    "binary_search",
]
