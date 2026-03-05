from __future__ import annotations

from pathlib import Path

import pytest

import SpliceGrapher.SpliceGraph as splice_graph_module
from SpliceGrapher.SpliceGraph import (
    GENE_REC,
    SpliceGraph,
    SpliceGraphNode,
    SpliceGraphParser,
    getFirstGraph,
    overlap,
)


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


def test_union_rejects_conflicting_known_strands_with_value_error() -> None:
    left = SpliceGraph("left", "chr1", "+")
    left.addNode("L1", 10, 20)

    right = SpliceGraph("right", "chr1", "-")
    right.addNode("R1", 30, 40)

    with pytest.raises(ValueError, match="conflicting strands"):
        left.union(right)


def test_union_requires_explicit_known_kwargs() -> None:
    left = SpliceGraph("left", "chr1", "+")
    left.addNode("L1", 10, 20)

    right = SpliceGraph("right", "chr1", "+")
    right.addNode("R1", 30, 40)

    merged = left.union(right, keepName=True)
    assert merged.getName() == "left"

    with pytest.raises(TypeError):
        left.union(right, keep_name=True)  # type: ignore[call-arg]


def test_splice_graph_node_defaults_do_not_share_parent_or_child_lists() -> None:
    left = SpliceGraphNode("L", 10, 20, "+", "chr1")
    right = SpliceGraphNode("R", 30, 40, "+", "chr1")
    parent = SpliceGraphNode("P", 1, 5, "+", "chr1")
    child = SpliceGraphNode("C", 50, 60, "+", "chr1")

    left.parents.append(parent)
    left.children.append(child)

    assert right.parents == []
    assert right.children == []


def _write_graph_header(path: Path, attrs: str) -> None:
    path.write_text(
        f"chr1\tSpliceGrapher\t{GENE_REC}\t1\t20\t.\t+\t.\t{attrs}\n",
        encoding="utf-8",
    )


def test_parser_preserves_equals_in_attribute_values(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1;Note=alpha=beta")

    parser = SpliceGraphParser(str(graph_path))
    graph = next(iter(parser.graphDict.values()))

    assert graph.attrs["Note"] == "alpha=beta"


def test_get_first_graph_uses_explicit_options(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1")

    graph = getFirstGraph(str(graph_path), annotate=False, verbose=False)
    assert graph.getName() == "G1"

    with pytest.raises(TypeError):
        getFirstGraph(str(graph_path), bad_option=True)  # type: ignore[call-arg]


def test_parser_progress_indicator_finishes_on_parse_error(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    graph_path = tmp_path / "invalid.gff3"
    _write_graph_header(graph_path, "ID")

    class DummyIndicator:
        finished = False

        def __init__(self, *_args, **_kwargs) -> None:
            pass

        def update(self) -> None:
            pass

        def finish(self) -> None:
            type(self).finished = True

    monkeypatch.setattr(splice_graph_module, "ProgressIndicator", DummyIndicator)

    with pytest.raises(ValueError, match="Illegal attribute field"):
        SpliceGraphParser(str(graph_path))

    assert DummyIndicator.finished is True


def test_parser_next_does_not_mask_internal_type_errors(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1")

    parser = SpliceGraphParser(str(graph_path))
    parser.graphId = "broken"  # type: ignore[assignment]

    with pytest.raises(TypeError):
        next(parser)


def test_parser_requires_explicit_constructor_kwargs(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1")

    with pytest.raises(TypeError):
        SpliceGraphParser(str(graph_path), bad_option=True)  # type: ignore[call-arg]


def test_overlap_helper_preserves_strict_boundary_behavior() -> None:
    graph = SpliceGraph("g", "chr1", "+")
    left = graph.addNode("L", 10, 20)
    touching = graph.addNode("T", 20, 30)
    crossing = graph.addNode("C", 19, 25)

    assert not overlap(left, touching)
    assert overlap(left, crossing)
