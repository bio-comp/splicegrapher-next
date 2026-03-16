"""Pure functions for annotating alternative splicing events on splice graphs."""

from __future__ import annotations

from typing import cast

from SpliceGrapher.core.enums import AlternativeSplicingEvent
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)
from SpliceGrapher.core.splice_graph.constants import AS_KEY
from SpliceGrapher.core.splice_graph.graph import SpliceGraph
from SpliceGrapher.core.splice_graph.node import SpliceGraphNode

EdgeBounds = tuple[int, int]


def _inner_edge_bounds(parent: SpliceGraphNode, child: SpliceGraphNode) -> EdgeBounds:
    positions = sorted([parent.minpos, parent.maxpos, child.minpos, child.maxpos])
    return positions[1], positions[2]


def _get_edge_bounds(graph: SpliceGraph) -> list[EdgeBounds]:
    """Return legacy inner edge bounds for all edges in the graph."""
    bounds: list[EdgeBounds] = []
    for parent_id, child_id in graph._nx_graph.edges:
        parent = cast(SpliceGraphNode, graph._nx_graph.nodes[parent_id]["data"])
        child = cast(SpliceGraphNode, graph._nx_graph.nodes[child_id]["data"])
        bounds.append(_inner_edge_bounds(parent, child))
    return bounds


def _annotate_alt3_events(nodes: list[SpliceGraphNode], graph: SpliceGraph) -> None:
    nodes_with_parents = [node for node in nodes if graph.predecessor_ids(node.id)]
    if not nodes_with_parents:
        return

    overlap_index = InMemoryIntervalIndex(nodes_with_parents)
    for node in nodes_with_parents:
        for other in overlap_index.overlaps(node, inclusive=False):
            if other == node or not intervals_overlap(node, other, inclusive=False):
                continue
            if node.acceptor_end() == other.acceptor_end():
                continue
            node.add_alt_form(AlternativeSplicingEvent.ALT3)
            other.add_alt_form(AlternativeSplicingEvent.ALT3)


def _annotate_alt5_events(nodes: list[SpliceGraphNode], graph: SpliceGraph) -> None:
    nodes_with_children = [node for node in nodes if graph.successor_ids(node.id)]
    if not nodes_with_children:
        return

    overlap_index = InMemoryIntervalIndex(nodes_with_children)
    for node in nodes_with_children:
        for other in overlap_index.overlaps(node, inclusive=False):
            if other == node or not intervals_overlap(node, other, inclusive=False):
                continue
            if node.donor_end() == other.donor_end():
                continue
            node.add_alt_form(AlternativeSplicingEvent.ALT5)
            other.add_alt_form(AlternativeSplicingEvent.ALT5)


def _annotate_branching_events(graph: SpliceGraph, nodes: list[SpliceGraphNode]) -> None:
    """Detect base ALT3/ALT5 motifs using graph-owned overlap queries."""
    _annotate_alt3_events(nodes, graph)
    _annotate_alt5_events(nodes, graph)


def _get_alt_site_event_groups(
    graph: SpliceGraph,
    nodes: list[SpliceGraphNode],
    event_type: AlternativeSplicingEvent,
) -> list[set[SpliceGraphNode]]:
    """Group overlapping ALT3 or ALT5 nodes using legacy containment exclusions."""
    event_nodes = sorted(
        [node for node in nodes if event_type in node.alt_form_set],
        key=lambda node: (node.minpos, node.maxpos, node.id),
    )
    if not event_nodes:
        return []

    def get_neighbors(node: SpliceGraphNode) -> set[SpliceGraphNode]:
        if event_type == AlternativeSplicingEvent.ALT3:
            return set(graph.predecessors(node))
        return set(graph.successors(node))

    overlap_index = InMemoryIntervalIndex(event_nodes)
    event_list: list[set[SpliceGraphNode]] = []
    stored: set[SpliceGraphNode] = set()

    for node in event_nodes:
        current_set = {node}
        neighbors = get_neighbors(node)

        for other in overlap_index.overlaps(node, inclusive=False):
            if other == node or not intervals_overlap(other, node, inclusive=False):
                continue

            other_neighbors = get_neighbors(other)
            node_contains_other_neighbors = (
                all(interval_contains(node, neighbor, strict=True) for neighbor in other_neighbors)
                if other_neighbors
                else False
            )
            other_contains_node_neighbors = (
                all(interval_contains(other, neighbor, strict=True) for neighbor in neighbors)
                if neighbors
                else False
            )

            if not (node_contains_other_neighbors or other_contains_node_neighbors):
                current_set.add(other)

        merged = False
        for existing_set in event_list:
            if existing_set & current_set:
                existing_set.update(current_set)
                stored.update(current_set)
                merged = True
                break

        if not merged and node not in stored:
            event_list.append(current_set)
            stored.update(current_set)

    return event_list


def _upgrade_to_alt_both(graph: SpliceGraph, nodes: list[SpliceGraphNode]) -> None:
    """Upgrade base ALT3/ALT5 events to ALTB3/ALTB5 when paths are unique."""
    alt5_groups = _get_alt_site_event_groups(graph, nodes, AlternativeSplicingEvent.ALT5)
    for group in alt5_groups:
        for node in group:
            node_acceptors = {child.acceptor_end() for child in graph.successors(node)}
            other_acceptors = {
                child.acceptor_end()
                for other in group
                if other != node
                for child in graph.successors(other)
            }
            if node_acceptors and not (node_acceptors & other_acceptors):
                node.remove_alt_form(AlternativeSplicingEvent.ALT5)
                node.add_alt_form(AlternativeSplicingEvent.ALTB5)

    alt3_groups = _get_alt_site_event_groups(graph, nodes, AlternativeSplicingEvent.ALT3)
    for group in alt3_groups:
        for node in group:
            node_donors = {parent.donor_end() for parent in graph.predecessors(node)}
            other_donors = {
                parent.donor_end()
                for other in group
                if other != node
                for parent in graph.predecessors(other)
            }
            if node_donors and not (node_donors & other_donors):
                node.remove_alt_form(AlternativeSplicingEvent.ALT3)
                node.add_alt_form(AlternativeSplicingEvent.ALTB3)


def annotate_graph_events(graph: SpliceGraph) -> None:
    """Annotate motif-level ES/IR/ALT3/ALT5 events on a networkx-backed graph."""
    nodes = graph.resolved_nodes()
    edge_bounds = _get_edge_bounds(graph)

    for node in nodes:
        node.alt_form_set.clear()
        node.attrs.pop(AS_KEY, None)

    for node in nodes:
        is_skipped = False
        is_retained = False
        for edge_min, edge_max in edge_bounds:
            if edge_min < node.minpos and node.maxpos < edge_max:
                is_skipped = True
            if node.minpos <= edge_min and edge_max <= node.maxpos:
                is_retained = True

        if is_skipped:
            if not graph.predecessor_ids(node.id):
                node.add_alt_form(AlternativeSplicingEvent.ALTI)
            elif not graph.successor_ids(node.id):
                node.add_alt_form(AlternativeSplicingEvent.ALTT)
            else:
                node.add_alt_form(AlternativeSplicingEvent.ES)

        if is_retained:
            node.add_alt_form(AlternativeSplicingEvent.IR)

    _annotate_branching_events(graph, nodes)
    _upgrade_to_alt_both(graph, nodes)


__all__ = ["annotate_graph_events"]
