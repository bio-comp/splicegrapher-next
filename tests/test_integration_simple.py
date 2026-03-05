from pathlib import Path

from SpliceGrapher.formats.annotation_io import load_gene_models
from SpliceGrapher.SpliceGraph import SpliceGraph, SpliceGraphParser
from tests.helpers.idiffir_fixture_builder import build_fixture


def _build_graph(name: str, start: int) -> SpliceGraph:
    graph = SpliceGraph(name, "chr1", "+")
    left = graph.addNode(f"{name}_left", start, start + 20)
    right = graph.addNode(f"{name}_right", start + 40, start + 60)
    left.addChild(right)
    return graph


def _load_single_graph(path: Path) -> SpliceGraph:
    parser = SpliceGraphParser(str(path))
    assert len(parser.graphDict) == 1
    return next(iter(parser.graphDict.values()))


def test_happy_path_load_union_write_cycle(tmp_path: Path) -> None:
    fixture = build_fixture(tmp_path)

    # Annotation path: load normalized gene model and ensure write roundtrip works.
    model = load_gene_models(str(fixture.gff3))
    model_out = tmp_path / "annotation.roundtrip.gff3"
    model.write_gff(str(model_out))
    model_text = model_out.read_text(encoding="utf-8")
    assert "\tgene\t" in model_text
    assert "\texon\t" in model_text

    # Graph path: load two serialized graphs, union them, then write merged output.
    left_path = tmp_path / "left.graph.gff3"
    right_path = tmp_path / "right.graph.gff3"
    _build_graph("left", 10).writeGFF(str(left_path))
    _build_graph("right", 15).writeGFF(str(right_path))

    left_graph = _load_single_graph(left_path)
    right_graph = _load_single_graph(right_path)

    merged = left_graph.union(right_graph)
    merged_out = tmp_path / "merged.graph.gff3"
    assert merged.writeGFF(str(merged_out))

    merged_loaded = _load_single_graph(merged_out)
    assert merged_loaded.getName()
    assert len(merged_loaded.nodeDict) >= 2
