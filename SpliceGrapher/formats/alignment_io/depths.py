"""Gene-depth helpers for alignment I/O."""

from __future__ import annotations

import numpy
import pysam
import structlog

from .collect import _record_strand
from .types import GeneBounds

LOGGER = structlog.get_logger(__name__)


def calculate_gene_depths(
    alignment_file: pysam.AlignmentFile,
    chromosome: str,
    gene: GeneBounds,
    *,
    margin: int = 0,
    verbose: bool = False,
) -> tuple[int, numpy.ndarray]:
    """Returns a relative start position and an array of read depths for a gene."""
    lo_bound = max(0, gene.minpos - margin)
    up_bound = gene.maxpos + margin

    spliced_reads = 0
    ungapped_reads = 0
    result = numpy.zeros(up_bound - lo_bound + 1, dtype=numpy.int32)

    for read in alignment_file.fetch(chromosome, lo_bound, up_bound):
        if read.is_unmapped:
            continue

        blocks = read.get_blocks()
        if len(blocks) > 1:
            if _record_strand(read) != gene.strand:
                continue
            spliced_reads += 1
        else:
            ungapped_reads += 1

        for block_start, block_end in blocks:
            start = max(block_start + 1, lo_bound) - lo_bound
            end = min(block_end, up_bound) - lo_bound
            if end > start:
                result[start:end] += 1

    if verbose:
        LOGGER.info(
            "alignment_depths_loaded",
            ungapped_reads=ungapped_reads,
            spliced_reads=spliced_reads,
            gene_id=gene.id,
        )

    return lo_bound, result
