"""Writer boundary for ``SpliceGraph`` serialization."""

from __future__ import annotations

from collections.abc import Iterable
from typing import TYPE_CHECKING, TextIO

import structlog

from SpliceGrapher.core.splice_graph import CHILD_REC, GENE_REC, PARENT_REC, SOURCE_NAME
from SpliceGrapher.shared.file_utils import open_output

if TYPE_CHECKING:
    from SpliceGrapher.core.splice_graph import SpliceGraph, SpliceGraphNode

LOGGER = structlog.get_logger(__name__)


def _gff_line_for_graph(graph: SpliceGraph) -> str:
    attrs_str = ";".join(f"{key}={value}" for key, value in sorted(graph.attrs.items()))
    return (
        f"{graph.chromosome}\t{SOURCE_NAME}\t{GENE_REC}\t"
        f"{graph.minpos}\t{graph.maxpos}\t.\t{graph.strand}\t.\t{attrs_str}"
    )


def _gff_line_for_node(
    node: SpliceGraphNode,
    parent_ids: Iterable[str] | None = None,
) -> str:
    parents = list(parent_ids) if parent_ids is not None else []
    record_type = CHILD_REC if parents else PARENT_REC
    attr_parts = [f"ID={node.id}"]
    if parents:
        attr_parts.append(f"Parent={','.join(parents)}")
    node_attrs = node.attribute_string()
    if node_attrs:
        attr_parts.append(node_attrs)
    attrs_str = ";".join(attr_parts)
    return (
        f"{node.chromosome}\t{SOURCE_NAME}\t{record_type}\t"
        f"{node.minpos}\t{node.maxpos}\t.\t{node.strand}\t.\t{attrs_str}"
    )


def write_splice_graph_gff(
    graph: SpliceGraph,
    out_path: str | TextIO,
    *,
    halt_on_error: bool = False,
) -> bool:
    reason = graph.validate()
    if reason:
        if halt_on_error:
            raise ValueError(
                f'Cannot write invalid splice graph {graph.get_name()} to file:\n"{reason}"\n'
            )
        LOGGER.warning(
            "writing_invalid_splice_graph",
            graph_name=graph.get_name(),
            reason=reason,
        )

    root_ids = {node_id for node_id, degree in graph._nx_graph.in_degree() if degree == 0}
    other_ids = set(graph._nx_graph.nodes) - root_ids

    def _node_sort_key(node_id: str) -> tuple[int, int, str]:
        node = graph._nx_graph.nodes[node_id]["data"]
        return (node.minpos, node.maxpos, node.id)

    with open_output(out_path) as out_stream:
        out_stream.write(f"{_gff_line_for_graph(graph)}\n")
        for node_id in sorted(root_ids, key=_node_sort_key):
            node = graph._nx_graph.nodes[node_id]["data"]
            out_stream.write(f"{_gff_line_for_node(node)}\n")
        for node_id in sorted(other_ids, key=_node_sort_key):
            node = graph._nx_graph.nodes[node_id]["data"]
            out_stream.write(f"{_gff_line_for_node(node, graph.predecessor_ids(node_id))}\n")
    return True


__all__ = ["write_splice_graph_gff"]
