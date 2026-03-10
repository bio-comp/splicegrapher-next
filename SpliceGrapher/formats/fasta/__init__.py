"""Modern FASTA reader and operation boundary."""

from __future__ import annotations

from .operations import (
    FastaRandomizer,
    fasta_count,
    fasta_split,
    get_sequence,
    truncate_sequences,
)
from .readers import FastaIterator, FastaSlice, fasta_get_by_name
from .records import FastaRecord, MalformedInput

__all__ = [
    "FastaIterator",
    "FastaRandomizer",
    "FastaRecord",
    "FastaSlice",
    "MalformedInput",
    "fasta_count",
    "fasta_get_by_name",
    "fasta_split",
    "get_sequence",
    "truncate_sequences",
]
