"""Alignment I/O public package boundary."""

from __future__ import annotations

from .api import (
    calculate_gene_depths,
    collect_alignment_data,
    read_alignment_chromosome_info,
    read_alignment_depths,
    read_alignment_headers,
    read_alignment_junctions,
    read_alignment_sequences,
    read_alignment_spans,
)
from .types import (
    AlignmentMap,
    AlignmentSource,
    ChromosomeInput,
    CollectResult,
    CollectResultWithAlignments,
    DepthMap,
    GeneBounds,
    JunctionMap,
    ReadDataSource,
    ReferencePath,
)

__all__ = [
    "AlignmentMap",
    "AlignmentSource",
    "ChromosomeInput",
    "CollectResult",
    "CollectResultWithAlignments",
    "DepthMap",
    "GeneBounds",
    "JunctionMap",
    "ReadDataSource",
    "ReferencePath",
    "calculate_gene_depths",
    "collect_alignment_data",
    "read_alignment_chromosome_info",
    "read_alignment_depths",
    "read_alignment_headers",
    "read_alignment_junctions",
    "read_alignment_sequences",
    "read_alignment_spans",
]
