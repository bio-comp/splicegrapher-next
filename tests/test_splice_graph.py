from __future__ import annotations

from pathlib import Path

import networkx as nx
import pytest

from SpliceGrapher.core.enums import (
    AlternativeSplicingEvent,
    AlternativeSplicingEventName,
    Strand,
)
from SpliceGrapher.core.splice_graph import GENE_REC, SpliceGraph, SpliceGraphNode
from SpliceGrapher.formats.parsers.splice_graph import SpliceGraphParser
from SpliceGrapher.formats.writers.splice_graph import write_splice_graph_gff


def _write_graph_header(path: Path, attrs: str) -> None:
    path.write_text(
        f"chr1\tSpliceGrapher\t{GENE_REC}\t1\t20\t.\t+\t.\t{attrs}\n",
        encoding="utf-8",
    )


def _load_single_graph(path: Path) -> SpliceGraph:
    parser = SpliceGraphParser(str(path))
    assert len(parser.graph_dict) == 1
    return next(iter(parser.graph_dict.values()))


def test_splice_graph_node_is_pure_data() -> None:
    node = SpliceGraphNode(
        id="exon_1",
        start=100,
        end=200,
        strand=Strand.PLUS,
        chrom="chr1",
    )

    assert node.id == "exon_1"
    assert node.minpos == 100
    assert node.maxpos == 200

    with pytest.raises(AttributeError):
        _ = node.parents  # type: ignore[attr-defined]

    with pytest.raises(AttributeError):
        _ = node.children  # type: ignore[attr-defined]

    with pytest.raises(AttributeError):
        _ = node.addChild  # type: ignore[attr-defined]

    with pytest.raises(AttributeError):
        _ = node.removeParent  # type: ignore[attr-defined]


def test_splice_graph_manages_topology() -> None:
    graph = SpliceGraph(name="test_gene", chromosome="chr1", strand=Strand.PLUS)

    node_a = graph.add_node(new_id="exon_1", start=100, end=200)
    node_b = graph.add_node(new_id="exon_2", start=300, end=400)

    assert isinstance(graph._nx_graph, nx.DiGraph)
    assert "exon_1" in graph._nx_graph.nodes
    assert graph._nx_graph.nodes["exon_1"]["data"] == node_a
    assert graph.node_dict["exon_2"] == node_b

    graph.add_edge("exon_1", "exon_2")

    assert graph._nx_graph.has_edge("exon_1", "exon_2")
    assert [node.id for node in graph.get_roots()] == ["exon_1"]
    assert [node.id for node in graph.get_leaves()] == ["exon_2"]


def test_splice_graph_missing_edge_handling() -> None:
    graph = SpliceGraph(name="test_gene", chromosome="chr1", strand=Strand.PLUS)
    graph.add_node(new_id="exon_1", start=100, end=200)

    with pytest.raises(ValueError, match="not found in graph"):
        graph.add_edge("exon_1", "phantom_exon")


def test_alt_form_set_normalizes_event_name_to_event_code() -> None:
    node = SpliceGraphNode("n1", 10, 20, "+", "chr1")

    node.add_alt_form(AlternativeSplicingEventName.ALT3.value)

    assert set(node.alt_forms()) == {AlternativeSplicingEvent.ALT3.value}


def test_parser_preserves_equals_in_attribute_values(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1;Note=alpha=beta")

    parser = SpliceGraphParser(str(graph_path))
    graph = next(iter(parser.graph_dict.values()))

    assert graph.attrs["Note"] == "alpha=beta"


def test_parser_iterates_over_graph_values_and_renames_loader(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1")

    parser = SpliceGraphParser(str(graph_path))

    assert not hasattr(parser, "loadFromFile")
    assert callable(parser.load_from_file)
    assert [graph.get_name() for graph in parser] == ["G1"]


def test_write_splice_graph_gff_roundtrip(tmp_path: Path) -> None:
    graph = SpliceGraph("graph_1", "chr1", Strand.PLUS)
    graph.add_node("exon_1", 100, 200)
    graph.add_node("exon_2", 300, 400)
    graph.add_edge("exon_1", "exon_2")

    with pytest.raises(AttributeError):
        _ = graph.writeGFF  # type: ignore[attr-defined]

    out_path = tmp_path / "graph.gff3"
    assert write_splice_graph_gff(graph, str(out_path)) is True

    loaded = _load_single_graph(out_path)

    assert set(loaded.node_dict) == {"exon_1", "exon_2"}
    assert [node.id for node in loaded.get_roots()] == ["exon_1"]
    assert [node.id for node in loaded.get_leaves()] == ["exon_2"]
