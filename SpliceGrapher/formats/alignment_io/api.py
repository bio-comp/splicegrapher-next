"""Public alignment I/O API over the internal alignment_io package modules."""

from __future__ import annotations

import sys
from typing import Literal, cast, overload

from SpliceGrapher.formats.depth_io import DepthSource

from . import collect as collect_ops
from .depths import calculate_gene_depths
from .sources import AlignmentStreamer
from .types import (
    AlignmentMap,
    AlignmentSource,
    ChromosomeInput,
    CollectResult,
    CollectResultWithAlignments,
    DepthMap,
    JunctionMap,
    ReadDataSource,
    ReferencePath,
)


def read_alignment_spans(
    source: AlignmentSource,
    *,
    chromosomes: ChromosomeInput = None,
    maxpos: int = sys.maxsize,
    reference_fasta: ReferencePath = None,
) -> AlignmentMap:
    """Read alignments and return ``(start, span)`` tuples per chromosome."""
    _, _, alignments = collect_ops._collect_pysam_data(
        source,
        alignments=True,
        chromosomes=chromosomes,
        junctions=False,
        maxpos=maxpos,
        reference_fasta=reference_fasta,
    )
    return alignments


def read_alignment_depths(
    source: ReadDataSource,
    *,
    verbose: bool = False,
    maxpos: int = sys.maxsize,
    chromosomes: ChromosomeInput = None,
    reference_fasta: ReferencePath = None,
) -> DepthMap:
    """Return read depths indexed by chromosome for SAM/BAM/CRAM sources."""
    if collect_ops._is_depths_source(source):
        depths, _ = collect_ops._collect_depths_source_data(
            cast(DepthSource, source),
            alignments=False,
            include_depths=True,
            maxpos=maxpos,
            junctions=False,
            minanchor=0,
            minjct=1,
            verbose=verbose,
        )
        return depths

    depths, _ = collect_ops._collect_pysam_data(
        cast(AlignmentSource, source),
        chromosomes=chromosomes,
        junctions=False,
        maxpos=maxpos,
        reference_fasta=reference_fasta,
    )
    return depths


def read_alignment_headers(
    source: AlignmentSource,
    *,
    reference_fasta: ReferencePath = None,
) -> list[str]:
    """Read a SAM/BAM/CRAM input and return the header strings."""
    streamer = AlignmentStreamer(source, reference_fasta=reference_fasta)
    with streamer.open_alignment() as sam_stream:
        return [line.strip() for line in str(sam_stream.header).splitlines() if line.strip()]


def read_alignment_chromosome_info(
    source: AlignmentSource,
    *,
    verbose: bool = False,
    reference_fasta: ReferencePath = None,
) -> set[str]:
    """Return chromosome names advertised by the alignment headers."""
    streamer = AlignmentStreamer(source, reference_fasta=reference_fasta)
    with streamer.open_alignment() as stream:
        chroms = {name for name in stream.references}

    if verbose:
        import structlog

        logger = structlog.get_logger(__name__)
        logger.info("sam_header_chromosome_count", chromosome_count=len(chroms))
        logger.info("sam_header_chromosomes", chromosomes=",".join(sorted(chroms)))

    return chroms


def read_alignment_junctions(
    source: ReadDataSource,
    *,
    verbose: bool = False,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    chromosomes: ChromosomeInput = None,
    reference_fasta: ReferencePath = None,
) -> JunctionMap:
    """Return splice junctions indexed by chromosome."""
    if collect_ops._is_depths_source(source):
        _, jcts = collect_ops._collect_depths_source_data(
            cast(DepthSource, source),
            alignments=False,
            include_depths=False,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            junctions=True,
            verbose=verbose,
        )
        return jcts

    _, junctions = collect_ops._collect_pysam_data(
        cast(AlignmentSource, source),
        chromosomes=chromosomes,
        maxpos=maxpos,
        minanchor=minanchor,
        minjct=minjct,
        reference_fasta=reference_fasta,
    )
    return junctions


@overload
def collect_alignment_data(
    source: ReadDataSource,
    *,
    include_alignments: Literal[True],
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    verbose: bool = False,
    reference_fasta: ReferencePath = None,
) -> CollectResultWithAlignments: ...


@overload
def collect_alignment_data(
    source: ReadDataSource,
    *,
    include_alignments: Literal[False] = False,
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    verbose: bool = False,
    reference_fasta: ReferencePath = None,
) -> CollectResult: ...


def collect_alignment_data(
    source: ReadDataSource,
    *,
    include_alignments: bool = False,
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    verbose: bool = False,
    reference_fasta: ReferencePath = None,
) -> CollectResult | CollectResultWithAlignments:
    """Return depths and junctions (and optionally alignments) from alignment input."""
    if collect_ops._is_depths_source(source):
        if include_alignments:
            return collect_ops._collect_depths_source_data(
                cast(DepthSource, source),
                alignments=True,
                include_depths=True,
                maxpos=maxpos,
                minanchor=minanchor,
                minjct=minjct,
                junctions=junctions,
                verbose=verbose,
            )
        return collect_ops._collect_depths_source_data(
            cast(DepthSource, source),
            alignments=False,
            include_depths=True,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            junctions=junctions,
            verbose=verbose,
        )

    _ = verbose
    if include_alignments:
        return collect_ops._collect_pysam_data(
            cast(AlignmentSource, source),
            alignments=True,
            chromosomes=chromosomes,
            junctions=junctions,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            reference_fasta=reference_fasta,
        )
    return collect_ops._collect_pysam_data(
        cast(AlignmentSource, source),
        alignments=False,
        chromosomes=chromosomes,
        junctions=junctions,
        maxpos=maxpos,
        minanchor=minanchor,
        minjct=minjct,
        reference_fasta=reference_fasta,
    )


def read_alignment_sequences(
    source: AlignmentSource,
    *,
    reference_fasta: ReferencePath = None,
) -> dict[str, str]:
    """Read a SAM/BAM/CRAM input and return sequence names and lengths."""
    streamer = AlignmentStreamer(source, reference_fasta=reference_fasta)
    with streamer.open_alignment() as stream:
        return {name: str(stream.get_reference_length(name)) for name in stream.references}


__all__ = [
    "calculate_gene_depths",
    "collect_alignment_data",
    "read_alignment_chromosome_info",
    "read_alignment_depths",
    "read_alignment_headers",
    "read_alignment_junctions",
    "read_alignment_sequences",
    "read_alignment_spans",
]
