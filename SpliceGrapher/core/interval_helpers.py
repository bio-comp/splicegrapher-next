"""Reusable interval overlap helpers and in-memory indexing utilities."""

from __future__ import annotations

import bisect
from collections.abc import Iterable, Sequence
from typing import Protocol, TypeVar, overload, runtime_checkable


@runtime_checkable
class IntervalBounds(Protocol):
    """Protocol for objects exposing interval start/end bounds."""

    @property
    def minpos(self) -> int: ...

    @property
    def maxpos(self) -> int: ...


_IntervalT = TypeVar("_IntervalT", bound=IntervalBounds)


def intervals_overlap(
    left: IntervalBounds,
    right: IntervalBounds,
    *,
    inclusive: bool = True,
) -> bool:
    """Return True when two interval bounds overlap."""
    if inclusive:
        return left.minpos <= right.maxpos and right.minpos <= left.maxpos
    return left.minpos < right.maxpos and right.minpos < left.maxpos


def interval_contains(
    container: IntervalBounds,
    candidate: IntervalBounds,
    *,
    strict: bool = False,
) -> bool:
    """Return True when ``container`` fully contains ``candidate``."""
    if strict:
        return container.minpos < candidate.minpos and container.maxpos > candidate.maxpos
    return container.minpos <= candidate.minpos and container.maxpos >= candidate.maxpos


class InMemoryIntervalIndex(Sequence[_IntervalT]):
    """Bisect-backed in-memory interval index over a sorted interval list."""

    def __init__(self, intervals: Sequence[_IntervalT]):
        if not intervals:
            raise ValueError("Cannot index an empty interval list")

        self._intervals = tuple(intervals)
        self._starts = tuple(interval.minpos for interval in self._intervals)

    @overload
    def __getitem__(self, index: int) -> _IntervalT: ...

    @overload
    def __getitem__(self, index: slice) -> Sequence[_IntervalT]: ...

    def __getitem__(self, index: int | slice) -> _IntervalT | Sequence[_IntervalT]:
        return self._intervals[index]

    def __len__(self) -> int:
        return len(self._intervals)

    def predecessor_or_containing(
        self,
        query: IntervalBounds,
        *,
        lo: int = 0,
        hi: int | None = None,
        overlap_window: int = 8,
    ) -> _IntervalT:
        """Return an interval that contains ``query`` or its immediate predecessor."""
        if hi is None:
            hi = len(self._intervals) - 1

        lo = max(0, lo)
        hi = min(hi, len(self._intervals) - 1)
        if lo > hi:
            raise ValueError("Invalid search bounds")

        idx = bisect.bisect_left(self._starts, query.minpos, lo=lo, hi=hi + 1)

        left = max(lo, idx - overlap_window)
        right = min(hi, idx + overlap_window)
        for i in range(left, right + 1):
            if interval_contains(self._intervals[i], query):
                return self._intervals[i]

        pred_idx = idx - 1
        if pred_idx < lo:
            return self._intervals[lo]
        if pred_idx > hi:
            return self._intervals[hi]
        return self._intervals[pred_idx]

    def overlaps(
        self,
        query: IntervalBounds,
        *,
        inclusive: bool = True,
    ) -> list[_IntervalT]:
        """Return all indexed intervals that overlap ``query``."""
        idx = bisect.bisect_left(self._starts, query.minpos)
        scan_start = max(0, idx - 1)
        result: list[_IntervalT] = []
        for interval in self._intervals[scan_start:]:
            if interval.minpos > query.maxpos:
                break
            if intervals_overlap(interval, query, inclusive=inclusive):
                result.append(interval)
        return result


def batch_overlaps(
    index: InMemoryIntervalIndex[_IntervalT],
    queries: Iterable[IntervalBounds],
    *,
    inclusive: bool = True,
) -> list[list[_IntervalT]]:
    """Bulk-overlap hook over a shared in-memory index."""
    return [index.overlaps(query, inclusive=inclusive) for query in queries]


__all__ = [
    "InMemoryIntervalIndex",
    "IntervalBounds",
    "batch_overlaps",
    "interval_contains",
    "intervals_overlap",
]
