"""Immutable genomic interval value object for gene-model domains."""

from __future__ import annotations

from dataclasses import dataclass

from .constants import SPLICE_DIMER_OFFSET, VALID_STRANDS


@dataclass(slots=True, frozen=True, eq=True)
class Locus:
    """Immutable genomic coordinate span."""

    chromosome: str
    strand: str
    minpos: int
    maxpos: int

    @classmethod
    def create(cls, chrom: str, start: int, end: int, strand: str) -> Locus:
        chrom_value = chrom.lower()
        strand_value = strand
        if strand_value not in VALID_STRANDS:
            raise ValueError(f"Unknown strand '{strand}'")
        return cls(
            chromosome=chrom_value,
            strand=strand_value,
            minpos=min(start, end),
            maxpos=max(start, end),
        )

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def contains(self, pos: int, strand: str) -> bool:
        return self.strand == strand and self.minpos <= pos <= self.maxpos

    def overlaps(self, other: Locus) -> bool:
        return (
            self.chromosome == other.chromosome
            and self.strand == other.strand
            and self.minpos <= other.maxpos
            and other.minpos <= self.maxpos
        )

    def start(self) -> int:
        return self.minpos if self.strand == "+" else self.maxpos

    def end(self) -> int:
        return self.maxpos if self.strand == "+" else self.minpos

    def acceptor(self, dimer_offset: int = SPLICE_DIMER_OFFSET) -> int:
        start = self.start()
        return start - dimer_offset if self.strand == "+" else start

    def donor(self, dimer_offset: int = SPLICE_DIMER_OFFSET) -> int:
        end = self.end()
        return end if self.strand == "+" else end - dimer_offset
