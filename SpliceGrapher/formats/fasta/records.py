"""FASTA record types and format-specific exceptions."""

from __future__ import annotations

from dataclasses import dataclass


class MalformedInput(Exception):
    """Raised when the input does not look like FASTA."""


@dataclass(frozen=True, slots=True)
class FastaRecord:
    """A single FASTA record."""

    header: str
    sequence: str

    def __str__(self) -> str:
        return f">{self.header}\n{self.sequence}\n"
