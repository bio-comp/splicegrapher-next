"""Pysam-backed collection helpers for alignment I/O."""

from __future__ import annotations

import os
import sys
from collections.abc import Mapping
from os import PathLike
from typing import Literal, overload

import numpy
import pysam

from SpliceGrapher.formats.depth_io import DepthSource, is_depths_file, read_depths
from SpliceGrapher.formats.junction import SpliceJunction, parse_junction_record

from .sources import AlignmentStreamer, _make_chromosome_set
from .types import (
    AlignmentMap,
    AlignmentSource,
    ChromosomeInput,
    CollectResult,
    CollectResultWithAlignments,
    DepthMap,
    JunctionKey,
    JunctionMap,
    ReadDataSource,
    ReferencePath,
)

STRAND_TAG = "XS"
JCT_CODE_TAG = "YC"
MATCH_CIGAR_OPS = frozenset({pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF})
QUERY_ONLY_CIGAR_OPS = frozenset({pysam.CINS, pysam.CSOFT_CLIP, pysam.CHARD_CLIP, pysam.CPAD})


def _next_match_anchor(cigar_tuples: list[tuple[int, int]], start_idx: int) -> int:
    """Return the contiguous downstream match length used as junction anchor."""
    anchor = 0
    for op, size in cigar_tuples[start_idx:]:
        if op in MATCH_CIGAR_OPS:
            anchor += size
            continue
        if op in QUERY_ONLY_CIGAR_OPS:
            continue
        break
    return anchor


def _build_splice_junctions(
    chromosome: str,
    start_pos: int,
    cigar_tuples: list[tuple[int, int]],
    strand: str,
    jct_code: str = "",
) -> list[SpliceJunction]:
    """Build splice junctions directly from pysam CIGAR tuples."""
    if not cigar_tuples:
        return []

    junctions = []
    reference_pos = start_pos
    left_anchor = 0

    for idx, (op, size) in enumerate(cigar_tuples):
        if op in MATCH_CIGAR_OPS:
            left_anchor += size
            reference_pos += size
            continue

        if op in QUERY_ONLY_CIGAR_OPS:
            continue

        if op == pysam.CDEL:
            reference_pos += size
            left_anchor = 0
            continue

        if op != pysam.CREF_SKIP:
            continue

        donor = reference_pos - 1
        reference_pos += size
        acceptor = reference_pos
        right_anchor = _next_match_anchor(cigar_tuples, idx + 1)

        if left_anchor > 0 and right_anchor > 0:
            junctions.append(
                SpliceJunction(
                    chromosome,
                    donor,
                    acceptor,
                    [left_anchor, right_anchor],
                    jct_code,
                    strand,
                )
            )

        left_anchor = 0

    return junctions


def _record_strand(record: pysam.AlignedSegment) -> str:
    if record.has_tag(STRAND_TAG):
        return str(record.get_tag(STRAND_TAG))
    return "-" if record.is_reverse else "+"


def _initial_depth_capacity(required_size: int, maxpos: int, reference_length: int | None) -> int:
    if maxpos < sys.maxsize:
        return maxpos + 1
    if reference_length is not None and reference_length >= 0:
        return max(required_size, reference_length + 2)
    return max(required_size, 1024)


def _ensure_depth_capacity(depth_array: numpy.ndarray, required_size: int) -> numpy.ndarray:
    if required_size <= depth_array.size:
        return depth_array

    new_size = max(required_size, depth_array.size * 2)
    expanded = numpy.zeros(new_size, dtype=depth_array.dtype)
    expanded[: depth_array.size] = depth_array
    return expanded


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
    """Bridge legacy depths-file inputs into the public alignment I/O contract."""
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


@overload
def _collect_pysam_data(
    source: AlignmentSource,
    *,
    alignments: Literal[True],
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    reference_fasta: ReferencePath = None,
) -> CollectResultWithAlignments: ...


@overload
def _collect_pysam_data(
    source: AlignmentSource,
    *,
    alignments: Literal[False] = False,
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    reference_fasta: ReferencePath = None,
) -> CollectResult: ...


def _collect_pysam_data(
    source: AlignmentSource,
    *,
    alignments: bool = False,
    chromosomes: ChromosomeInput = None,
    junctions: bool = True,
    maxpos: int = sys.maxsize,
    minanchor: int = 0,
    minjct: int = 1,
    reference_fasta: ReferencePath = None,
) -> CollectResult | CollectResultWithAlignments:
    """Collect depths/junctions/alignment spans using pysam-backed iteration."""
    chrom_set = _make_chromosome_set(chromosomes)

    align: AlignmentMap = {}
    depths: DepthMap = {}
    junction_tmp: dict[JunctionKey, SpliceJunction] = {}
    limit: dict[str, int] = {}

    streamer = AlignmentStreamer(source, reference_fasta=reference_fasta)

    with streamer.open_alignment() as alignment:
        for rec in alignment.fetch(until_eof=True):
            if rec.is_unmapped:
                continue

            chrom_name = alignment.get_reference_name(rec.reference_id)
            if chrom_name is None:
                continue

            chrom = chrom_name.lower()
            if chrom_set and chrom not in chrom_set:
                continue

            start_pos = rec.reference_start + 1
            if start_pos > maxpos:
                continue

            reference_end = (
                rec.reference_end if rec.reference_end is not None else rec.reference_start
            )
            required_size = reference_end + 2

            if chrom not in depths:
                reference_length: int | None = None
                try:
                    reference_length = alignment.get_reference_length(chrom_name)
                except (KeyError, ValueError):
                    reference_length = None
                depths[chrom] = numpy.zeros(
                    _initial_depth_capacity(required_size, maxpos, reference_length),
                    dtype=numpy.int32,
                )
                limit[chrom] = 0
                if alignments:
                    align[chrom] = []
            elif maxpos == sys.maxsize:
                depths[chrom] = _ensure_depth_capacity(depths[chrom], required_size)

            if alignments and chrom not in align:
                align[chrom] = []

            for block_start, block_end in rec.get_blocks():
                block_lo = block_start + 1
                block_hi = block_end + 1

                if block_lo > maxpos:
                    break

                if maxpos == sys.maxsize:
                    depths[chrom] = _ensure_depth_capacity(depths[chrom], block_hi)
                    slice_hi = block_hi
                else:
                    slice_hi = min(block_hi, maxpos + 1)

                if slice_hi <= block_lo:
                    continue

                depths[chrom][block_lo:slice_hi] += 1

                if alignments:
                    align[chrom].append((block_lo, slice_hi - block_lo))

            limit[chrom] = max(limit[chrom], reference_end + 2)

            if not junctions or not rec.cigartuples:
                continue

            strand = _record_strand(rec)
            jct_code = str(rec.get_tag(JCT_CODE_TAG)) if rec.has_tag(JCT_CODE_TAG) else ""
            jct_list = _build_splice_junctions(chrom, start_pos, rec.cigartuples, strand, jct_code)
            if not jct_list:
                continue

            for new_jct in jct_list:
                key = (chrom, strand, new_jct.p1, new_jct.p2)
                existing = junction_tmp.get(key)
                if existing is None:
                    junction_tmp[key] = new_jct
                else:
                    existing.update(new_jct)

    normalized_depths: DepthMap = {}
    for chrom, depth_array in depths.items():
        if maxpos < sys.maxsize:
            keep = min(maxpos + 1, depth_array.size)
        else:
            keep = min(limit[chrom], depth_array.size)
        normalized_depths[chrom] = depth_array[:keep].copy()

    normalized_junctions: JunctionMap = {}
    if junctions:
        for (chrom, _strand, _anchor, _acceptor), jct in sorted(junction_tmp.items()):
            if jct.count >= minjct and jct.min_anchor() >= minanchor:
                normalized_junctions.setdefault(chrom, []).append(jct)

    if alignments:
        return normalized_depths, normalized_junctions, align
    return normalized_depths, normalized_junctions


__all__ = [
    "_collect_depths_source_data",
    "_collect_pysam_data",
    "_depth_map_to_arrays",
    "_is_depths_source",
    "_record_strand",
]
