from __future__ import annotations

from dataclasses import dataclass

from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    batch_overlaps,
    interval_contains,
    intervals_overlap,
)


@dataclass(slots=True)
class _Interval:
    minpos: int
    maxpos: int


def test_intervals_overlap_supports_inclusive_and_exclusive_modes() -> None:
    left = _Interval(10, 20)
    right = _Interval(20, 30)

    assert intervals_overlap(left, right, inclusive=True)
    assert not intervals_overlap(left, right, inclusive=False)


def test_interval_contains_supports_strict_and_non_strict_modes() -> None:
    container = _Interval(10, 20)
    strict_child = _Interval(11, 19)
    boundary_child = _Interval(10, 20)

    assert interval_contains(container, strict_child, strict=True)
    assert not interval_contains(container, boundary_child, strict=True)
    assert interval_contains(container, boundary_child, strict=False)


def test_in_memory_interval_index_returns_containing_or_predecessor() -> None:
    intervals = [_Interval(1, 10), _Interval(20, 30), _Interval(40, 50)]
    index = InMemoryIntervalIndex(intervals)

    contained = index.predecessor_or_containing(_Interval(21, 22))
    predecessor = index.predecessor_or_containing(_Interval(32, 33))

    assert contained == intervals[1]
    assert predecessor == intervals[1]


def test_batch_overlaps_uses_shared_interval_index() -> None:
    intervals = [_Interval(1, 10), _Interval(20, 30), _Interval(40, 50)]
    index = InMemoryIntervalIndex(intervals)
    queries = [_Interval(2, 4), _Interval(25, 45)]

    result = batch_overlaps(index, queries, inclusive=True)

    assert result == [[intervals[0]], [intervals[1], intervals[2]]]
