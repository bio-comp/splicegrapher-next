"""Alignment I/O helpers for SAM/BAM/CRAM inputs."""

import io
import os
import re
import sys
from collections.abc import Iterator, Mapping
from contextlib import contextmanager
from os import PathLike
from typing import Literal, Protocol, TypeAlias, cast, overload

import numpy
import pysam
import structlog

from SpliceGrapher.formats.depth_io import DepthSource, is_depths_file, read_depths
from SpliceGrapher.formats.junction import SpliceJunction, parse_junction_record

LOGGER = structlog.get_logger(__name__)

STRAND_TAG = "XS"
JCT_CODE_TAG = "YC"
NULL_CHROMOSOME = "*"

CIGAR_TOKEN = re.compile(r"[0-9]+[A-Z]")

MATCH_CIGAR_OPS = frozenset({pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF})
QUERY_ONLY_CIGAR_OPS = frozenset({pysam.CINS, pysam.CSOFT_CLIP, pysam.CHARD_CLIP, pysam.CPAD})

ChromosomeInput: TypeAlias = str | list[str] | tuple[str, ...] | set[str] | None
ReferencePath: TypeAlias = str | PathLike[str] | None
SamLine: TypeAlias = str | bytes
SamLineCollection: TypeAlias = list[SamLine] | tuple[SamLine, ...] | set[SamLine]
SamTextSource: TypeAlias = io.IOBase | SamLineCollection
AlignmentPath: TypeAlias = str | PathLike[str]
AlignmentSource: TypeAlias = AlignmentPath | SamTextSource
ReadDataSource: TypeAlias = AlignmentSource | DepthSource
DepthMap: TypeAlias = dict[str, numpy.ndarray]
JunctionMap: TypeAlias = dict[str, list[SpliceJunction]]
AlignmentMap: TypeAlias = dict[str, list[tuple[int, int]]]
JunctionKey: TypeAlias = tuple[str, str, int, int]
CollectResult: TypeAlias = tuple[DepthMap, JunctionMap]
CollectResultWithAlignments: TypeAlias = tuple[DepthMap, JunctionMap, AlignmentMap]


class GeneBounds(Protocol):
    id: str
    strand: str
    minpos: int
    maxpos: int


def _is_bam_path(file_path: str | PathLike[str]) -> bool:
    """Return ``True`` when the given path looks like a BAM file."""
    return str(file_path).lower().endswith(".bam")


def _is_cram_path(file_path: str | PathLike[str]) -> bool:
    """Return ``True`` when the given path looks like a CRAM file."""
    return str(file_path).lower().endswith(".cram")


def _make_chromosome_set(chromosomes: ChromosomeInput) -> set[str] | None:
    """Convert a chromosome selector into a lowercase set."""
    if not chromosomes:
        return None
    if isinstance(chromosomes, str):
        chromosomes = [chromosomes]
    return {chrom.lower() for chrom in chromosomes}


# Source normalization and alignment opening helpers.
def _is_alignment_path(source: AlignmentSource) -> bool:
    """Return ``True`` when ``source`` is a filesystem path to an alignment file."""
    return isinstance(source, (str, PathLike)) and os.path.isfile(os.fspath(source))


def _as_text_line(line: str | bytes) -> str:
    return line.decode("utf-8") if isinstance(line, bytes) else line


def _iter_sam_lines(source: SamTextSource) -> Iterator[str]:
    if isinstance(source, io.IOBase):
        for line in source:
            yield _as_text_line(line)
        return

    if isinstance(source, (list, tuple, set)):
        for sam_line in source:
            yield _as_text_line(sam_line)
        return

    raise ValueError(f"Unrecognized SAM input source: {type(source)!r}")


def _reference_consumed_from_cigar(cigar: str) -> int:
    total = 0
    for token in CIGAR_TOKEN.findall(cigar):
        if token[-1] in {"M", "D", "N", "=", "X"}:
            total += int(token[:-1])
    return total


def _synthesize_sq_headers(records: list[str]) -> list[str]:
    lengths: dict[str, int] = {}
    for record in records:
        fields = record.rstrip("\n").split("\t")
        if len(fields) < 6:
            continue

        chromosome = fields[2]
        if chromosome == NULL_CHROMOSOME:
            continue

        try:
            start = int(fields[3])
        except ValueError:
            continue

        ref_len = _reference_consumed_from_cigar(fields[5])
        if ref_len <= 0:
            seq_len = len(fields[9]) if len(fields) > 9 else 0
            ref_len = seq_len

        end = start + max(ref_len - 1, 0)
        lengths[chromosome] = max(lengths.get(chromosome, 0), end)

    return [f"@SQ\tSN:{chrom}\tLN:{length}\n" for chrom, length in sorted(lengths.items())]


@contextmanager
def _open_alignment_source(
    source: AlignmentSource,
    *,
    reference_fasta: ReferencePath = None,
) -> Iterator[pysam.AlignmentFile]:
    """Open paths directly and wrap in-memory SAM lines in a BytesIO buffer."""
    if _is_alignment_path(source):
        with _open_alignment_file(str(source), reference_fasta=reference_fasta) as alignment:
            yield alignment
        return

    if isinstance(source, (str, PathLike)):
        raise ValueError(f"Alignment source path not found: {source}")

    headers: list[str] = []
    records: list[str] = []
    for raw_line in _iter_sam_lines(source):
        line = raw_line if raw_line.endswith("\n") else f"{raw_line}\n"
        if line.startswith("@"):
            headers.append(line)
        elif line.strip():
            records.append(line)

    if not any(line.startswith("@SQ") for line in headers):
        headers.extend(_synthesize_sq_headers(records))
    if not any(line.startswith("@HD") for line in headers):
        headers.insert(0, "@HD\tVN:1.6\tSO:unknown\n")

    memory_buffer = io.BytesIO()
    try:
        for header in headers:
            memory_buffer.write(header.encode("utf-8"))
        for record in records:
            memory_buffer.write(record.encode("utf-8"))
        memory_buffer.seek(0)
        with pysam.AlignmentFile(memory_buffer, "r") as alignment:
            yield alignment
    finally:
        memory_buffer.close()


def _open_alignment_file(
    path: str | PathLike[str],
    *,
    reference_fasta: ReferencePath = None,
) -> pysam.AlignmentFile:
    """Open SAM/BAM/CRAM input with pysam and CRAM reference safeguards."""
    path_string = os.fspath(path)

    if _is_bam_path(path_string):
        return pysam.AlignmentFile(path_string, "rb")

    if _is_cram_path(path_string):
        try:
            if reference_fasta:
                return pysam.AlignmentFile(
                    path_string,
                    "rc",
                    reference_filename=str(reference_fasta),
                )
            return pysam.AlignmentFile(path_string, "rc")
        except ValueError as exc:
            if "reference" in str(exc).lower():
                raise ValueError(
                    "Unable to decode CRAM without reference. "
                    "Re-run with reference_fasta=<path-to-fasta>."
                ) from exc
            raise

    return pysam.AlignmentFile(path_string, "r")


# Pysam-backed collection helpers.
class AlignmentStreamer:
    """Unified alignment streamer for SAM/BAM/CRAM and in-memory SAM sources."""

    def __init__(
        self,
        source: AlignmentSource,
        *,
        reference_fasta: ReferencePath = None,
    ) -> None:
        self.source = source
        self.reference_fasta = reference_fasta

    @contextmanager
    def open_alignment(self) -> Iterator[pysam.AlignmentFile]:
        with _open_alignment_source(self.source, reference_fasta=self.reference_fasta) as alignment:
            yield alignment


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


def _depth_map_to_arrays(
    depths: Mapping[str, list[int] | numpy.ndarray],
) -> dict[str, numpy.ndarray]:
    """Normalize depth maps to int32 NumPy arrays."""
    return {
        chrom: numpy.asarray(chrom_depths, dtype=numpy.int32)
        for chrom, chrom_depths in depths.items()
    }


# Depths-file fallback bridge.
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
            if jct.count >= minjct and jct.minAnchor() >= minanchor:
                normalized_junctions.setdefault(chrom, []).append(jct)

    if alignments:
        return normalized_depths, normalized_junctions, align

    return normalized_depths, normalized_junctions


#
# Public alignment API.
#
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


def read_alignment_spans(
    source: AlignmentSource,
    *,
    chromosomes: ChromosomeInput = None,
    maxpos: int = sys.maxsize,
    reference_fasta: ReferencePath = None,
) -> AlignmentMap:
    """Read alignments and return ``(start, span)`` tuples per chromosome."""
    _, _, alignments = _collect_pysam_data(
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

    depths, _ = _collect_pysam_data(
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
    """Reads a SAM file and returns just the header strings as a list."""
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
        LOGGER.info("sam_header_chromosome_count", chromosome_count=len(chroms))
        LOGGER.info("sam_header_chromosomes", chromosomes=",".join(sorted(chroms)))

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

    _, junctions = _collect_pysam_data(
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
        return _collect_pysam_data(
            cast(AlignmentSource, source),
            alignments=True,
            chromosomes=chromosomes,
            junctions=junctions,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            reference_fasta=reference_fasta,
        )
    return _collect_pysam_data(
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
    """Reads a SAM/BAM/CRAM input and returns sequence names and lengths."""
    streamer = AlignmentStreamer(source, reference_fasta=reference_fasta)
    with streamer.open_alignment() as stream:
        return {name: str(stream.get_reference_length(name)) for name in stream.references}


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
