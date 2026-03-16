"""Graph model for the splice-graph core."""

from __future__ import annotations

import sys
from collections.abc import Iterable

import networkx as nx
import structlog

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import Strand

from .constants import END_CODON_KEY, ID_ATTR, SOURCE_NAME, START_CODON_KEY
from .node import NullNode, SpliceGraphNode

LOGGER = structlog.get_logger(__name__)


class SpliceGraph:
    """NetworkX-backed manager for splice graph topology."""

    def __init__(self, name: str, chromosome: str, strand: str | Strand) -> None:
        self.chromosome = chromosome
        self.strand = coerce_enum(strand, Strand, field="strand").value
        self.attrs: dict[str, str] = {}
        self.minpos = sys.maxsize
        self.maxpos = 0
        self._nx_graph: nx.DiGraph = nx.DiGraph()
        self.set_name(name)

    @property
    def node_dict(self) -> dict[str, SpliceGraphNode]:
        return {
            node_id: node_data["data"] for node_id, node_data in self._nx_graph.nodes(data=True)
        }

    def _sorted_nodes(self, node_ids: Iterable[str]) -> list[SpliceGraphNode]:
        nodes = [self._nx_graph.nodes[node_id]["data"] for node_id in node_ids]
        return sorted(nodes, key=lambda node: (node.minpos, node.maxpos, node.id))

    def _find_existing_node(self, start: int, end: int) -> SpliceGraphNode | None:
        probe = NullNode(start, end)
        for node in self.node_dict.values():
            if node == probe:
                return node
        return None

    def _recompute_bounds(self) -> None:
        nodes = list(self.node_dict.values())
        if not nodes:
            self.minpos = sys.maxsize
            self.maxpos = 0
            return
        self.minpos = min(node.minpos for node in nodes)
        self.maxpos = max(node.maxpos for node in nodes)

    def _node(self, node_or_id: str | SpliceGraphNode) -> SpliceGraphNode:
        if isinstance(node_or_id, SpliceGraphNode):
            return node_or_id
        return self._nx_graph.nodes[node_or_id]["data"]

    def add_node(self, new_id: str, start: int, end: int) -> SpliceGraphNode:
        existing = self._find_existing_node(start, end)
        if existing is not None:
            return existing
        node = SpliceGraphNode(new_id, start, end, self.strand, self.chromosome)
        self._nx_graph.add_node(node.id, data=node)
        self.minpos = min(self.minpos, node.minpos)
        self.maxpos = max(self.maxpos, node.maxpos)
        return node

    def add_codons(self, codon_list: Iterable[tuple[int, int]], codon_type: str) -> None:
        for codon in codon_list:
            for node in self.node_dict.values():
                node.add_codon(codon, codon_type)

    def add_edge(self, pid: str, cid: str) -> None:
        if pid not in self._nx_graph:
            raise ValueError(f"Error adding edge: parent node {pid} not found in graph")
        if cid not in self._nx_graph:
            raise ValueError(f"Error adding edge: child node {cid} not found in graph")
        self._nx_graph.add_edge(pid, cid)

    def add_end_codons(self, codon_list: Iterable[tuple[int, int]]) -> None:
        self.add_codons(codon_list, END_CODON_KEY)

    def add_start_codons(self, codon_list: Iterable[tuple[int, int]]) -> None:
        self.add_codons(codon_list, START_CODON_KEY)

    def adjust(self, adjustment: int) -> None:
        for node in self.node_dict.values():
            node.update(node.minpos + adjustment, node.maxpos + adjustment)
        if self.node_dict:
            self.minpos += adjustment
            self.maxpos += adjustment

    def attribute_string(self) -> str:
        return ";".join(f"{key}={self.attrs[key]}" for key in sorted(self.attrs))

    def delete_node(self, node_id: str | SpliceGraphNode) -> SpliceGraphNode:
        node = self._node(node_id)
        self._nx_graph.remove_node(node.id)
        self._recompute_bounds()
        return node

    def get_leaves(self) -> list[SpliceGraphNode]:
        return self._sorted_nodes(
            node_id for node_id, degree in self._nx_graph.out_degree() if degree == 0
        )

    def get_name(self) -> str:
        return self.attrs[ID_ATTR]

    def get_node(self, start: int, end: int) -> SpliceGraphNode | None:
        return self._find_existing_node(start, end)

    def get_roots(self) -> list[SpliceGraphNode]:
        return self._sorted_nodes(
            node_id for node_id, degree in self._nx_graph.in_degree() if degree == 0
        )

    def is_empty(self) -> bool:
        return self._nx_graph.number_of_nodes() == 0

    def predecessors(self, node_or_id: str | SpliceGraphNode) -> list[SpliceGraphNode]:
        node = self._node(node_or_id)
        return self._sorted_nodes(self._nx_graph.predecessors(node.id))

    def predecessor_ids(self, node_id: str) -> list[str]:
        return sorted(self._nx_graph.predecessors(node_id))

    def resolved_nodes(self) -> list[SpliceGraphNode]:
        return [node for node in self.node_dict.values() if not node.is_unresolved()]

    def set_name(self, name: str) -> None:
        self.attrs[ID_ATTR] = name

    def successors(self, node_or_id: str | SpliceGraphNode) -> list[SpliceGraphNode]:
        node = self._node(node_or_id)
        return self._sorted_nodes(self._nx_graph.successors(node.id))

    def successor_ids(self, node_id: str) -> list[str]:
        return sorted(self._nx_graph.successors(node_id))

    def unresolved_nodes(self) -> list[SpliceGraphNode]:
        return [node for node in self.node_dict.values() if node.is_unresolved()]

    def validate(self, halt: bool = False) -> str | None:
        reason: str | None = None
        if self._nx_graph.number_of_nodes() == 0:
            reason = "Graph is empty."
        else:
            roots = self.get_roots()
            leaves = self.get_leaves()
            if not roots or not leaves:
                reason = (
                    f"Graph is missing roots or leaves ({len(roots)} roots, {len(leaves)} leaves)"
                )
            elif not nx.is_directed_acyclic_graph(self._nx_graph):
                reason = "Graph contains cycles (not a valid DAG)."

        if reason and halt:
            LOGGER.error("illegal_graph_validation_failed", reason=reason)
            raise ValueError(f"Illegal graph:\n{reason}")
        return reason

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1 if self.maxpos >= self.minpos else 0

    def __str__(self) -> str:
        return (
            f"{self.get_name()} ({self.strand}) {self.minpos}-{self.maxpos} "
            f"({self._nx_graph.number_of_nodes()} nodes)"
        )


__all__ = ["SOURCE_NAME", "SpliceGraph"]
