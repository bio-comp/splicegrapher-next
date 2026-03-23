"""Public alignment I/O API over the internal alignment_io package modules."""

from __future__ import annotations

import os
import sys
from collections.abc import Mapping
from os import PathLike
from typing import Literal, cast, overload

import numpy

from SpliceGrapher.formats.depth_io import DepthSource, is_depths_file, read_depths
from SpliceGrapher.formats.junction import parse_junction_record

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


def _depth_map_to_arrays(depths: Mapping[str, list[int] | numpy.ndarray]) -> DepthMap:
    """Normalize depth maps to int32 NumPy arrays."""
    return {
        chrom: numpy.asarray(chrom_depths, dtype=numpy.int32)
        for chrom, chrom_depths in depths.items()
    }


def _is_depths_source(source: ReadDataSource) -> bool:
    """Return ``True`` when ``source`` should be handled as a depths file."""
    if not isinstance(source, (str, PathLike)):
        return False

    source_path = os.fspath(source)
    lower = source_path.lower()
    if lower.endswith((".sam", ".bam", ".cram")):
        return False

    if os.path.isfile(source_path):
        try:
            with open(source_path, "rb") as stream:
                magic = stream.read(4)
        except OSError:
            magic = b""
        if magic in {b"BAM\x01", b"CRAM"}:
            return False

    return is_depths_file(source_path)


@overload
def _collect_depths_source_data(
    source: DepthSource,
    *,
    alignments: Literal[True],
    include_depths: bool = True,
    junctions: bool,
    maxpos: int,
    minanchor: int,
    minjct: int,
    verbose: bool,
) -> CollectResultWithAlignments: ...


@overload
def _collect_depths_source_data(
    source: DepthSource,
    *,
    alignments: Literal[False] = False,
    include_depths: bool = True,
    junctions: bool,
    maxpos: int,
    minanchor: int,
    minjct: int,
    verbose: bool,
) -> CollectResult: ...


def _collect_depths_source_data(
    source: DepthSource,
    *,
    alignments: bool = False,
    include_depths: bool = True,
    junctions: bool,
    maxpos: int,
    minanchor: int,
    minjct: int,
    verbose: bool,
) -> CollectResult | CollectResultWithAlignments:
    """Bridge depths-file inputs into the public alignment I/O contract."""
    depth_map, jcts = read_depths(
        source,
        parse_junction=parse_junction_record,
        maxpos=maxpos,
        minanchor=minanchor,
        minjct=minjct,
        depths=include_depths,
        junctions=junctions,
        verbose=verbose,
    )
    normalized_depths = _depth_map_to_arrays(depth_map)
    if alignments:
        empty_alignments: AlignmentMap = {}
        return normalized_depths, jcts, empty_alignments
    return normalized_depths, jcts


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
    if _is_depths_source(source):
        depths, _ = _collect_depths_source_data(
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
    if _is_depths_source(source):
        _, jcts = _collect_depths_source_data(
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
    if _is_depths_source(source):
        if include_alignments:
            return _collect_depths_source_data(
                cast(DepthSource, source),
                alignments=True,
                include_depths=True,
                maxpos=maxpos,
                minanchor=minanchor,
                minjct=minjct,
                junctions=junctions,
                verbose=verbose,
            )
        return _collect_depths_source_data(
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
