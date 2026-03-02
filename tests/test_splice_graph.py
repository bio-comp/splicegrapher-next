from __future__ import annotations

import pytest

from SpliceGrapher.SpliceGraph import SpliceGraph, overlap


def test_union_uses_resolved_strand_when_merging_unknown_strand_graph() -> None:
    left = SpliceGraph("left", "chr1", ".")
    left.addNode("L1", 10, 20)

    right = SpliceGraph("right", "chr1", "+")
    right.addNode("R1", 30, 40)

    merged = left.union(right)

    assert merged.strand == "+"


def test_union_rejects_invalid_runtime_strand_value() -> None:
    left = SpliceGraph("left", "chr1", "+")
    left.addNode("L1", 10, 20)

    right = SpliceGraph("right", "chr1", "+")
    right.addNode("R1", 30, 40)
    right.strand = "?"

    with pytest.raises(ValueError):
        left.union(right)


def test_overlap_helper_preserves_strict_boundary_behavior() -> None:
    graph = SpliceGraph("g", "chr1", "+")
    left = graph.addNode("L", 10, 20)
    touching = graph.addNode("T", 20, 30)
    crossing = graph.addNode("C", 19, 25)

    assert not overlap(left, touching)
    assert overlap(left, crossing)
