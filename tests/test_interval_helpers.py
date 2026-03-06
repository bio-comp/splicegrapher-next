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


def test_in_memory_interval_index_defensively_sorts_unsorted_inputs() -> None:
    intervals = [_Interval(40, 50), _Interval(1, 10), _Interval(20, 30)]
    index = InMemoryIntervalIndex(intervals)

    sorted_bounds = [(interval.minpos, interval.maxpos) for interval in index]
    overlaps = index.overlaps(_Interval(5, 25))

    assert sorted_bounds == [(1, 10), (20, 30), (40, 50)]
    assert [(interval.minpos, interval.maxpos) for interval in overlaps] == [(1, 10), (20, 30)]


def test_batch_overlaps_uses_shared_interval_index() -> None:
    intervals = [_Interval(1, 10), _Interval(20, 30), _Interval(40, 50)]
    index = InMemoryIntervalIndex(intervals)
    queries = [_Interval(2, 4), _Interval(25, 45)]

    result = batch_overlaps(index, queries, inclusive=True)

    assert result == [[intervals[0]], [intervals[1], intervals[2]]]


def test_in_memory_interval_index_overlap_parity_handles_nested_long_intervals() -> None:
    intervals = [
        _Interval(1, 100),
        _Interval(50, 60),
        _Interval(70, 80),
        _Interval(90, 95),
        _Interval(110, 120),
    ]
    query = _Interval(92, 93)
    index = InMemoryIntervalIndex(intervals)

    legacy = [
        interval
        for interval in intervals
        if interval.minpos <= query.maxpos and query.minpos <= interval.maxpos
    ]

    assert index.overlaps(query) == legacy
