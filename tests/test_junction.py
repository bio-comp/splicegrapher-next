from __future__ import annotations

import pytest

from SpliceGrapher.formats.junction import SpliceJunction, parse_junction_record


def test_junction_round_trips_via_depth_record_string() -> None:
    junction = SpliceJunction("chr1", 10, 20, (2, 3), "K", "+")
    junction.count = 7

    parsed = parse_junction_record(junction.to_string())

    assert parsed.chromosome == "chr1"
    assert parsed.minpos == 10
    assert parsed.maxpos == 20
    assert tuple(parsed.anchors) == (2, 3)
    assert parsed.sj_code == "K"
    assert parsed.count == 7
    assert parsed.strand == "+"


def test_junction_objects_sort_by_chromosome_and_coordinates() -> None:
    junctions = [
        SpliceJunction("chr2", 10, 20, (2, 2), "K", "+"),
        SpliceJunction("chr1", 30, 40, (2, 2), "K", "+"),
        SpliceJunction("chr1", 10, 20, (2, 2), "K", "+"),
    ]

    ordered = sorted(junctions)

    assert [j.chromosome for j in ordered] == ["chr1", "chr1", "chr2"]
    assert [j.minpos for j in ordered] == [10, 30, 10]


def test_junction_update_merges_counts_and_max_anchor_values() -> None:
    left = SpliceJunction("chr1", 10, 20, (2, 3), "K", "+")
    right = SpliceJunction("chr1", 10, 20, (4, 1), "K", "+")
    right.count = 5

    left.update(right)

    assert left.count == 6
    assert tuple(left.anchors) == (4, 3)


def test_junction_update_mismatch_raises_value_error() -> None:
    left = SpliceJunction("chr1", 10, 20, (2, 2), "K", "+")
    other = SpliceJunction("chr1", 12, 20, (2, 2), "K", "+")

    with pytest.raises(ValueError, match="Splice sites do not match"):
        left.update(other)
