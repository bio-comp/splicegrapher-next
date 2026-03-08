"""Pure functions for comparing and subtracting splice graphs."""

from __future__ import annotations

from typing import cast

from SpliceGrapher.core.splice_graph import SpliceGraph, SpliceGraphNode

EdgeCoords = tuple[int, int, int, int]


def _get_edge_coords(graph: SpliceGraph) -> set[EdgeCoords]:
    """Extract canonical edge coordinates for cross-graph comparisons."""
    coords: set[EdgeCoords] = set()
    for parent_id, child_id in graph._nx_graph.edges:
        parent = cast(SpliceGraphNode, graph._nx_graph.nodes[parent_id]["data"])
        child = cast(SpliceGraphNode, graph._nx_graph.nodes[child_id]["data"])
        coords.add((parent.minpos, parent.maxpos, child.minpos, child.maxpos))
    return coords


def jaccard_coefficients(graph_a: SpliceGraph, graph_b: SpliceGraph) -> tuple[float, float]:
    """Return Jaccard coefficients for resolved nodes and coordinate-based edges."""
    nodes_a = set(graph_a.resolvedNodes())
    nodes_b = set(graph_b.resolvedNodes())
    node_union = nodes_a | nodes_b
    node_coeff = (len(nodes_a & nodes_b) / len(node_union)) if node_union else 1.0

    edges_a = _get_edge_coords(graph_a)
    edges_b = _get_edge_coords(graph_b)
    edge_union = edges_a | edges_b
    edge_coeff = (len(edges_a & edges_b) / len(edge_union)) if edge_union else 1.0
    return float(node_coeff), float(edge_coeff)


def recall(graph_a: SpliceGraph, graph_b: SpliceGraph) -> tuple[float, float]:
    """Return recall values using ``graph_b`` as the denominator set."""
    nodes_b = set(graph_b.resolvedNodes())
    node_coeff = (len(set(graph_a.resolvedNodes()) & nodes_b) / len(nodes_b)) if nodes_b else 1.0

    edges_b = _get_edge_coords(graph_b)
    edge_coeff = (len(_get_edge_coords(graph_a) & edges_b) / len(edges_b)) if edges_b else 1.0
    return float(node_coeff), float(edge_coeff)


def equivalent_graphs(graph_a: SpliceGraph, graph_b: SpliceGraph) -> bool:
    """Return ``True`` when two graphs have identical node/edge coordinate sets."""
    if len(graph_a.nodeDict) != len(graph_b.nodeDict):
        return False
    return set(graph_a.nodeDict.values()) == set(graph_b.nodeDict.values()) and _get_edge_coords(
        graph_a
    ) == _get_edge_coords(graph_b)


__all__ = ["equivalent_graphs", "jaccard_coefficients", "recall"]
