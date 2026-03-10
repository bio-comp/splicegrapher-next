"""Higher-level FASTA operations."""

from __future__ import annotations

import random
import sys
from os import PathLike
from pathlib import Path
from typing import TextIO

import structlog

from SpliceGrapher.formats.fasta.readers import FastaIterator, FastaSource
from SpliceGrapher.formats.fasta.records import FastaRecord
from SpliceGrapher.shared.file_utils import open_output

LOGGER = structlog.get_logger(__name__)


class FastaRandomizer:
    """Load all FASTA sequences and return random subsets."""

    def __init__(self, fasta_file: FastaSource):
        self.records = list(FastaIterator(fasta_file))
        self.record_range = range(len(self.records))

    def random_records(self, count: int) -> list[FastaRecord]:
        if count > len(self.records):
            raise ValueError(f"{count} FASTA records requested; only {len(self.records)} stored.")
        record_indices = random.sample(self.record_range, count)
        return [self.records[index] for index in record_indices]

    def __str__(self) -> str:
        return f"FastaRandomizer - {len(self.records)} records"


def get_sequence(source: FastaSource, name: str) -> FastaRecord | None:
    """Return the record in ``source`` with the given header."""
    return FastaIterator(source)[name]


def fasta_count(source: FastaSource) -> int:
    """Count FASTA records in a supported source."""
    return sum(1 for _record in FastaIterator(source))


def fasta_split(
    file_name: str | PathLike[str],
    num_files: int,
    directory: str | PathLike[str] | None = None,
) -> None:
    """Split a FASTA file into the requested number of output files."""
    source_path = Path(file_name)
    num_records = fasta_count(source_path)

    if directory is None:
        base = source_path.with_suffix("")
        suffix = source_path.suffix
    else:
        base = Path(directory) / source_path.stem
        suffix = source_path.suffix

    record_number = 1
    file_number = 0
    records_per_file = float(num_records) / float(num_files)
    record_limit = 0.0
    output_handle: TextIO | None = None

    for record in FastaIterator(source_path):
        if record_number > round(record_limit):
            if output_handle is not None:
                output_handle.close()
            file_number += 1
            output_path = Path(f"{base}.{file_number}{suffix}")
            output_handle = output_path.open("w", encoding="utf-8")
            record_limit += records_per_file
        assert output_handle is not None
        output_handle.write(str(record))
        record_number += 1

    if output_handle is not None:
        output_handle.close()


def truncate_sequences(
    fasta_file: FastaSource,
    exon_size: int,
    intron_size: int,
    *,
    acceptor: bool = False,
    output_file: str | Path | TextIO | None = None,
    verbose: bool = False,
) -> None:
    """Truncate FASTA sequences around their midpoint."""
    if verbose:
        LOGGER.info("loading_fasta_sequences", source=str(fasta_file))

    output_target: str | Path | TextIO = sys.stdout if output_file is None else output_file
    if verbose and output_file is not None:
        LOGGER.info("writing_truncated_fasta", output=str(output_file))

    with open_output(output_target) as output_handle:
        for record in FastaIterator(fasta_file):
            midpoint = len(record.sequence) // 2
            if acceptor:
                truncated = FastaRecord(
                    record.header,
                    record.sequence[midpoint - intron_size : midpoint + exon_size],
                )
            else:
                truncated = FastaRecord(
                    record.header,
                    record.sequence[midpoint - exon_size : midpoint + intron_size],
                )
            output_handle.write(str(truncated))
