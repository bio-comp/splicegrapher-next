from __future__ import annotations

from SpliceGrapher.SpliceGraph import SpliceGraph


def test_union_uses_resolved_strand_when_merging_unknown_strand_graph() -> None:
    left = SpliceGraph("left", "chr1", ".")
    left.addNode("L1", 10, 20)

    right = SpliceGraph("right", "chr1", "+")
    right.addNode("R1", 30, 40)

    merged = left.union(right)

    assert merged.strand == "+"
