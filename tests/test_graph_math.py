from __future__ import annotations

from SpliceGrapher.core.graph_math import (
    equivalent_graphs,
    jaccard_coefficients,
    recall,
)
from SpliceGrapher.core.splice_graph import SpliceGraph


def _build_linear_graph(
    name: str,
    left_id: str,
    right_id: str,
    *,
    with_tail: bool = False,
) -> SpliceGraph:
    graph = SpliceGraph(name, "chr1", "+")
    graph.addNode(left_id, 10, 20)
    graph.addNode(right_id, 30, 40)
    graph.addEdge(left_id, right_id)
    if with_tail:
        graph.addNode(f"{name}_tail", 50, 60)
        graph.addEdge(right_id, f"{name}_tail")
    return graph


def test_equivalent_graphs_ignore_node_ids_and_graph_names() -> None:
    left = _build_linear_graph("left", "left_a", "left_b")
    right = _build_linear_graph("right", "right_x", "right_y")

    assert equivalent_graphs(left, right) is True


def test_jaccard_coefficients_compare_edges_by_coordinates() -> None:
    left = _build_linear_graph("left", "left_a", "left_b")
    right = _build_linear_graph("right", "right_x", "right_y", with_tail=True)

    node_coeff, edge_coeff = jaccard_coefficients(left, right)

    assert node_coeff == 2.0 / 3.0
    assert edge_coeff == 1.0 / 2.0


def test_recall_uses_graph_b_as_denominator_for_nodes_and_edges() -> None:
    left = _build_linear_graph("left", "left_a", "left_b")
    right = _build_linear_graph("right", "right_x", "right_y", with_tail=True)

    node_coeff, edge_coeff = recall(left, right)

    assert node_coeff == 2.0 / 3.0
    assert edge_coeff == 1.0 / 2.0
