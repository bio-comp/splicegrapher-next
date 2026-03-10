"""Alignment source normalization and opening helpers."""

from __future__ import annotations

import io
import os
import re
from collections.abc import Iterator
from contextlib import contextmanager
from os import PathLike
from typing import cast

import pysam

from .types import (
    AlignmentPath,
    AlignmentSource,
    ChromosomeInput,
    ReferencePath,
    SamTextSource,
)

NULL_CHROMOSOME = "*"
CIGAR_TOKEN = re.compile(r"[0-9]+[A-Z]")


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
        path_source = cast(AlignmentPath, source)
        with _open_alignment_file(path_source, reference_fasta=reference_fasta) as alignment:
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
