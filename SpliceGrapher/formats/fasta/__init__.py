"""Modern FASTA reader and operation boundary."""

from __future__ import annotations

from .records import FastaRecord, MalformedInput

__all__ = [
    "FastaRecord",
    "MalformedInput",
]
