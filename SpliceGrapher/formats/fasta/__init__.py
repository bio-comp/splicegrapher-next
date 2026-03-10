"""Modern FASTA reader and operation boundary."""

from __future__ import annotations

from .readers import FastaIterator, FastaSlice, fasta_get_by_name
from .records import FastaRecord, MalformedInput

__all__ = [
    "FastaIterator",
    "FastaRecord",
    "FastaSlice",
    "MalformedInput",
    "fasta_get_by_name",
]
