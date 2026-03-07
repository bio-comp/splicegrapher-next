from pathlib import Path

from SpliceGrapher.formats.annotation_io import load_gene_models
from SpliceGrapher.formats.writers.splice_graph import write_splice_graph_gff
from SpliceGrapher.SpliceGraph import SpliceGraph, SpliceGraphParser
from tests.helpers.idiffir_fixture_builder import build_fixture


def _build_graph(name: str, start: int) -> SpliceGraph:
    graph = SpliceGraph(name, "chr1", "+")
    graph.addNode(f"{name}_left", start, start + 20)
    graph.addNode(f"{name}_right", start + 40, start + 60)
    graph.addEdge(f"{name}_left", f"{name}_right")
    return graph


def _load_single_graph(path: Path) -> SpliceGraph:
    parser = SpliceGraphParser(str(path))
    assert len(parser.graphDict) == 1
    return next(iter(parser.graphDict.values()))


def test_happy_path_load_and_write_cycle(tmp_path: Path) -> None:
    fixture = build_fixture(tmp_path)

    model = load_gene_models(str(fixture.gff3))
    model_out = tmp_path / "annotation.roundtrip.gff3"
    model.write_gff(str(model_out))
    model_text = model_out.read_text(encoding="utf-8")
    assert "\tgene\t" in model_text
    assert "\texon\t" in model_text

    graph_path = tmp_path / "graph.gff3"
    graph = _build_graph("left", 10)
    assert write_splice_graph_gff(graph, str(graph_path)) is True

    loaded = _load_single_graph(graph_path)
    assert loaded.getName() == "left"
    assert set(loaded.nodeDict) == {"left_left", "left_right"}
    assert [node.id for node in loaded.getRoots()] == ["left_left"]
    assert [node.id for node in loaded.getLeaves()] == ["left_right"]
