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
    assert len(parser.graphDict) == 1
    return next(iter(parser.graphDict.values()))


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

    node_a = graph.addNode(newId="exon_1", start=100, end=200)
    node_b = graph.addNode(newId="exon_2", start=300, end=400)

    assert isinstance(graph._nx_graph, nx.DiGraph)
    assert "exon_1" in graph._nx_graph.nodes
    assert graph._nx_graph.nodes["exon_1"]["data"] == node_a
    assert graph.nodeDict["exon_2"] == node_b

    graph.addEdge("exon_1", "exon_2")

    assert graph._nx_graph.has_edge("exon_1", "exon_2")
    assert [node.id for node in graph.getRoots()] == ["exon_1"]
    assert [node.id for node in graph.getLeaves()] == ["exon_2"]


def test_splice_graph_missing_edge_handling() -> None:
    graph = SpliceGraph(name="test_gene", chromosome="chr1", strand=Strand.PLUS)
    graph.addNode(newId="exon_1", start=100, end=200)

    with pytest.raises(ValueError, match="not found in graph"):
        graph.addEdge("exon_1", "phantom_exon")


def test_alt_form_set_normalizes_event_name_to_event_code() -> None:
    node = SpliceGraphNode("n1", 10, 20, "+", "chr1")

    node.addAltForm(AlternativeSplicingEventName.ALT3.value)

    assert set(node.altForms()) == {AlternativeSplicingEvent.ALT3.value}


def test_parser_preserves_equals_in_attribute_values(tmp_path: Path) -> None:
    graph_path = tmp_path / "graph.gff3"
    _write_graph_header(graph_path, "ID=G1;Note=alpha=beta")

    parser = SpliceGraphParser(str(graph_path))
    graph = next(iter(parser.graphDict.values()))

    assert graph.attrs["Note"] == "alpha=beta"


def test_write_splice_graph_gff_roundtrip(tmp_path: Path) -> None:
    graph = SpliceGraph("graph_1", "chr1", Strand.PLUS)
    graph.addNode("exon_1", 100, 200)
    graph.addNode("exon_2", 300, 400)
    graph.addEdge("exon_1", "exon_2")

    with pytest.raises(AttributeError):
        _ = graph.writeGFF  # type: ignore[attr-defined]

    out_path = tmp_path / "graph.gff3"
    assert write_splice_graph_gff(graph, str(out_path)) is True

    loaded = _load_single_graph(out_path)

    assert set(loaded.nodeDict) == {"exon_1", "exon_2"}
    assert [node.id for node in loaded.getRoots()] == ["exon_1"]
    assert [node.id for node in loaded.getLeaves()] == ["exon_2"]
