from collections import defaultdict
from pathlib import Path

import numpy
import pysam

# TODO(#41): Remove this temporary parity shim after Phase 1 parser migration
# stabilizes. Removal criteria:
# 1) #36 merged; 2) CI green for two consecutive weeks (or one stable patch cycle);
# 3) fixture-based golden outputs replace code-based legacy oracle checks.


MATCH = 0
INSERT = 1
DELETE = 2
GAP = 3


def _open_alignment(path: Path, reference_fasta: Path | None = None) -> pysam.AlignmentFile:
    suffix = path.suffix.lower()
    if suffix == ".sam":
        return pysam.AlignmentFile(path, "r")
    if suffix == ".cram":
        if reference_fasta:
            return pysam.AlignmentFile(path, "rc", reference_filename=str(reference_fasta))
        return pysam.AlignmentFile(path, "rc")
    return pysam.AlignmentFile(path, "rb")


def _legacy_process_read(
    read: pysam.AlignedSegment,
    depths: numpy.ndarray,
    junctions: defaultdict[tuple[int, int], int],
    min_idx: int,
) -> None:
    pos = read.pos
    for op, length in read.cigar:
        if op == MATCH:
            start = max(0, pos - min_idx)
            end = min(len(depths), max(0, (pos + length) - min_idx))
            for idx in range(start, end):
                depths[idx] += 1
            pos += length
        elif op == INSERT:
            continue
        elif op == DELETE:
            # alignment_io does not increment depths across deletion spans.
            pos += length
        elif op == GAP:
            junctions[(pos - min_idx, (pos + length) - min_idx)] += 1
            pos += length


def compute_depths_and_junctions(
    alignment_path: Path,
    chrom: str,
    start: int,
    end: int,
    *,
    reference_fasta: Path | None = None,
) -> tuple[numpy.ndarray, dict[tuple[int, int], int]]:
    min_idx = start - 1
    depths = numpy.zeros(end - start + 1, int)
    junctions: defaultdict[tuple[int, int], int] = defaultdict(int)

    with _open_alignment(alignment_path, reference_fasta=reference_fasta) as handle:
        try:
            reads = handle.fetch(chrom, start, end)
            for read in reads:
                _legacy_process_read(read, depths, junctions, min_idx)
        except (ValueError, OSError):
            for read in handle:
                if read.reference_name != chrom:
                    continue
                if read.reference_end is None:
                    continue
                if read.reference_end < start or read.reference_start > end:
                    continue
                _legacy_process_read(read, depths, junctions, min_idx)

    return depths, dict(junctions)
