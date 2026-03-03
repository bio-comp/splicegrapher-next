"""Alignment I/O helpers for SAM/BAM/CRAM inputs."""

import io
import os
import re
import sys
from collections.abc import Iterator
from contextlib import contextmanager

import numpy
import pysam
import structlog

from SpliceGrapher.formats.shortread_compat import (
    SpliceJunction,
    is_depths_file,
    read_depths,
)

LOGGER = structlog.get_logger(__name__)

STRAND_TAG = "XS"
JCT_CODE_TAG = "YC"
NULL_CHROMOSOME = "*"

CIGAR_TOKEN = re.compile(r"[0-9]+[A-Z]")

MATCH_CIGAR_OPS = frozenset({pysam.CMATCH, pysam.CEQUAL, pysam.CDIFF})
QUERY_ONLY_CIGAR_OPS = frozenset({pysam.CINS, pysam.CSOFT_CLIP, pysam.CHARD_CLIP, pysam.CPAD})


def isBamFile(filePath):
    """Simple heuristic returns True if path is a BAM file; false otherwise."""
    return str(filePath).lower().endswith(".bam")


def isCramFile(filePath):
    """Simple heuristic returns True if path is a CRAM file; false otherwise."""
    return str(filePath).lower().endswith(".cram")


def makeChromosomeSet(chromList):
    """Converts a string or a list of chromosomes into a set of unique values."""
    chrom_set = None
    if chromList:
        if isinstance(chromList, str):
            chromList = [chromList]
        chrom_set = {chrom.lower() for chrom in chromList}
    return chrom_set


def _is_alignment_path(source):
    """Return ``True`` when ``source`` is a filesystem path to an alignment file."""
    return isinstance(source, str) and os.path.isfile(source)


def _as_text_line(line: str | bytes) -> str:
    return line.decode("utf-8") if isinstance(line, bytes) else line


def _iter_sam_lines(source: object) -> Iterator[str]:
    if isinstance(source, io.IOBase):
        for line in source:
            yield _as_text_line(line)
        return

    if isinstance(source, (list, tuple, set)):
        for line in source:
            if not isinstance(line, (str, bytes)):
                raise ValueError(f"Unsupported SAM line type: {type(line)!r}")
            yield _as_text_line(line)
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
def _open_alignment_source(source: object, *, reference_fasta=None):
    """Open paths directly and wrap in-memory SAM lines in a BytesIO buffer."""
    if _is_alignment_path(source):
        with _open_alignment_file(str(source), reference_fasta=reference_fasta) as alignment:
            yield alignment
        return

    if isinstance(source, str):
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


def _open_alignment_file(path, *, reference_fasta=None):
    """Open SAM/BAM/CRAM input with pysam and CRAM reference safeguards."""
    if isBamFile(path):
        return pysam.AlignmentFile(path, "rb")

    if isCramFile(path):
        try:
            if reference_fasta:
                return pysam.AlignmentFile(path, "rc", reference_filename=str(reference_fasta))
            return pysam.AlignmentFile(path, "rc")
        except ValueError as exc:
            if "reference" in str(exc).lower():
                raise ValueError(
                    "Unable to decode CRAM without reference. "
                    "Re-run with reference_fasta=<path-to-fasta>."
                ) from exc
            raise

    return pysam.AlignmentFile(path, "r")


def _is_depths_source(source) -> bool:
    """Return ``True`` when ``source`` should be handled as a depths file."""
    if not isinstance(source, str):
        return False

    lower = source.lower()
    if lower.endswith((".sam", ".bam", ".cram")):
        return False

    if os.path.isfile(source):
        try:
            with open(source, "rb") as stream:
                magic = stream.read(4)
        except OSError:
            magic = b""
        if magic in {b"BAM\x01", b"CRAM"}:
            return False

    return is_depths_file(source)


class AlignmentStreamer:
    """Unified alignment streamer for SAM/BAM/CRAM and in-memory SAM sources."""

    def __init__(self, source, *, reference_fasta=None):
        self.source = source
        self.reference_fasta = reference_fasta

    @contextmanager
    def open_alignment(self):
        with _open_alignment_source(self.source, reference_fasta=self.reference_fasta) as alignment:
            yield alignment

    def stream(self, *, chromosomes=None):
        chrom_set = makeChromosomeSet(chromosomes)
        with self.open_alignment() as alignment:
            for rec in alignment.fetch(until_eof=True):
                if rec.is_unmapped:
                    continue

                if chrom_set:
                    chrom_name = alignment.get_reference_name(rec.reference_id)
                    if not chrom_name or chrom_name.lower() not in chrom_set:
                        continue

                yield rec


def _next_match_anchor(cigar_tuples, start_idx):
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


def _build_splice_junctions(chromosome, start_pos, cigar_tuples, strand, jct_code=""):
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


def _record_strand(record) -> str:
    if record.has_tag(STRAND_TAG):
        return record.get_tag(STRAND_TAG)
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


def _depth_map_to_arrays(depths: dict[str, list[int] | numpy.ndarray]) -> dict[str, numpy.ndarray]:
    """Normalize depth maps to int32 NumPy arrays."""
    return {
        chrom: numpy.asarray(chrom_depths, dtype=numpy.int32)
        for chrom, chrom_depths in depths.items()
    }


def _collect_pysam_data(
    source,
    *,
    alignments=False,
    chromosomes=None,
    junctions=True,
    maxpos=sys.maxsize,
    minanchor=0,
    minjct=1,
    reference_fasta=None,
):
    """Collect depths/junctions/alignment spans using pysam-backed iteration."""
    chrom_set = makeChromosomeSet(chromosomes)

    align = {}
    depths: dict[str, numpy.ndarray] = {}
    junction_tmp = {}
    limit = {}

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
                reference_length = None
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
            jct_code = rec.get_tag(JCT_CODE_TAG) if rec.has_tag(JCT_CODE_TAG) else ""
            jct_list = _build_splice_junctions(chrom, start_pos, rec.cigartuples, strand, jct_code)
            if not jct_list:
                continue

            junction_tmp.setdefault(chrom, {})
            junction_tmp[chrom].setdefault(strand, {})
            for new_jct in jct_list:
                junction_tmp[chrom][strand].setdefault(new_jct.p1, {})
                junction_tmp[chrom][strand][new_jct.p1].setdefault(new_jct.p2, None)
                existing = junction_tmp[chrom][strand][new_jct.p1][new_jct.p2]
                if existing is None:
                    junction_tmp[chrom][strand][new_jct.p1][new_jct.p2] = new_jct
                else:
                    existing.update(new_jct)

    normalized_depths = {}
    for chrom, depth_array in depths.items():
        if maxpos < sys.maxsize:
            keep = min(maxpos + 1, depth_array.size)
        else:
            keep = min(limit[chrom], depth_array.size)

        normalized_depths[chrom] = depth_array[:keep].copy()

    normalized_junctions = {}
    if junctions:
        for chrom in sorted(junction_tmp.keys()):
            normalized_junctions[chrom] = []
            for strand in sorted(junction_tmp[chrom].keys()):
                for anchor in sorted(junction_tmp[chrom][strand].keys()):
                    for acceptor in sorted(junction_tmp[chrom][strand][anchor].keys()):
                        jct = junction_tmp[chrom][strand][anchor][acceptor]
                        if jct.count >= minjct and jct.minAnchor() >= minanchor:
                            normalized_junctions[chrom].append(jct)

    if alignments:
        return normalized_depths, normalized_junctions, align

    return normalized_depths, normalized_junctions


def pysamStrand(pysamRecord, tagDict=None):
    """Convenience method returns the strand given by the pysam record."""
    if tagDict and STRAND_TAG in tagDict:
        return tagDict[STRAND_TAG]
    return _record_strand(pysamRecord)


def pysamReadDepths(bamFile, chromosome, gene, *, margin=0, verbose=False):
    """Returns a relative start position and an array of read depths for a gene."""
    loBound = max(0, gene.minpos - margin)
    upBound = gene.maxpos + margin

    nSpliced = 0
    nUngapped = 0
    result = numpy.zeros(upBound - loBound + 1, dtype=numpy.int32)

    for read in bamFile.fetch(chromosome, loBound, upBound):
        if read.is_unmapped:
            continue

        blocks = read.get_blocks()
        if len(blocks) > 1:
            if pysamStrand(read) != gene.strand:
                continue
            nSpliced += 1
        else:
            nUngapped += 1

        for block_start, block_end in blocks:
            start = max(block_start + 1, loBound) - loBound
            end = min(block_end, upBound) - loBound
            if end > start:
                result[start:end] += 1

    if verbose:
        LOGGER.info(
            "pysam_read_depths_loaded",
            ungapped_reads=nUngapped,
            spliced_reads=nSpliced,
            gene_id=gene.id,
        )

    return loBound, result


def getSamAlignments(
    samRecords,
    *,
    verbose=False,
    chromosomes=None,
    maxpos=sys.maxsize,
    reference_fasta=None,
):
    """Read alignments and return ``(start, span)`` tuples per chromosome."""
    _ = verbose
    _, _, alignments = _collect_pysam_data(
        samRecords,
        alignments=True,
        chromosomes=chromosomes,
        junctions=False,
        maxpos=maxpos,
        reference_fasta=reference_fasta,
    )
    return alignments


def getSamDepths(
    samRecords,
    *,
    verbose=False,
    maxpos=sys.maxsize,
    chromosomes=None,
    reference_fasta=None,
):
    """Return read depths indexed by chromosome for SAM/BAM/CRAM sources."""
    if _is_depths_source(samRecords):
        depths, _ = read_depths(
            samRecords,
            maxpos=maxpos,
            junctions=False,
            verbose=verbose,
        )
        return _depth_map_to_arrays(depths)

    depths, _ = _collect_pysam_data(
        samRecords,
        chromosomes=chromosomes,
        junctions=False,
        maxpos=maxpos,
        reference_fasta=reference_fasta,
    )
    return depths


def getSamHeaders(samRecords, *, reference_fasta=None):
    """Reads a SAM file and returns just the header strings as a list."""
    streamer = AlignmentStreamer(samRecords, reference_fasta=reference_fasta)
    with streamer.open_alignment() as sam_stream:
        return [line.strip() for line in str(sam_stream.header).splitlines() if line.strip()]


def getSamHeaderInfo(samStream, *, verbose=False, reference_fasta=None):
    """Return chromosome names from headers and no seed alignment line."""
    streamer = AlignmentStreamer(samStream, reference_fasta=reference_fasta)
    with streamer.open_alignment() as stream:
        chroms = {name for name in stream.references}

    if verbose:
        LOGGER.info("sam_header_chromosome_count", chromosome_count=len(chroms))
        LOGGER.info("sam_header_chromosomes", chromosomes=",".join(sorted(chroms)))

    return chroms, None


def getSamJunctions(
    samRecords,
    *,
    verbose=False,
    maxpos=sys.maxsize,
    minanchor=0,
    minjct=1,
    chromosomes=None,
    reference_fasta=None,
):
    """Return splice junctions indexed by chromosome."""
    if _is_depths_source(samRecords):
        _, jcts = read_depths(
            samRecords,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            depths=False,
            verbose=verbose,
        )
        return jcts

    _, junctions = _collect_pysam_data(
        samRecords,
        chromosomes=chromosomes,
        maxpos=maxpos,
        minanchor=minanchor,
        minjct=minjct,
        reference_fasta=reference_fasta,
    )
    return junctions


def getSamReadData(
    samRecords,
    *,
    alignments=False,
    chromosomes=None,
    junctions=True,
    maxpos=sys.maxsize,
    minanchor=0,
    minjct=1,
    verbose=False,
    reference_fasta=None,
):
    """Return depths and junctions (and optionally alignments) from alignment input."""
    if _is_depths_source(samRecords):
        depths, jcts = read_depths(
            samRecords,
            maxpos=maxpos,
            minanchor=minanchor,
            minjct=minjct,
            depths=True,
            junctions=junctions,
            verbose=verbose,
        )
        return _depth_map_to_arrays(depths), jcts

    _ = verbose
    return _collect_pysam_data(
        samRecords,
        alignments=alignments,
        chromosomes=chromosomes,
        junctions=junctions,
        maxpos=maxpos,
        minanchor=minanchor,
        minjct=minjct,
        reference_fasta=reference_fasta,
    )


def getSamSequences(samRecords, *, reference_fasta=None):
    """Reads a SAM/BAM/CRAM input and returns sequence names and lengths."""
    streamer = AlignmentStreamer(samRecords, reference_fasta=reference_fasta)
    with streamer.open_alignment() as stream:
        return {name: str(stream.get_reference_length(name)) for name in stream.references}
