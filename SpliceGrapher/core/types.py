"""Typed interval value objects and protocols shared across SpliceGrapher."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Protocol, runtime_checkable

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import Strand


@runtime_checkable
class IntervalCoordinates(Protocol):
    """Protocol contract for overlap-capable coordinate objects."""

    @property
    def chromosome(self) -> str: ...

    @property
    def minpos(self) -> int: ...

    @property
    def maxpos(self) -> int: ...

    @property
    def strand(self) -> Strand | str: ...


@dataclass(frozen=True, slots=True, order=True)
class GenomicInterval:
    """Immutable, ordered interval value object."""

    chromosome: str
    start: int
    end: int
    strand: Strand = Strand.UNKNOWN

    def __post_init__(self) -> None:
        chromosome = self.chromosome.strip()
        if not chromosome:
            raise ValueError("chromosome must be non-empty")
        if self.start < 1 or self.end < 1:
            raise ValueError("interval coordinates must be >= 1")
        if self.start > self.end:
            raise ValueError("interval start must be <= end")

        object.__setattr__(self, "chromosome", chromosome)
        object.__setattr__(self, "strand", coerce_enum(self.strand, Strand, field="strand"))

    @classmethod
    def from_raw(
        cls,
        chromosome: str,
        start: int,
        end: int,
        strand: Strand | str = Strand.UNKNOWN,
    ) -> GenomicInterval:
        """Create an interval from raw parser-style values."""
        return cls(
            chromosome=chromosome,
            start=start,
            end=end,
            strand=coerce_enum(strand, Strand, field="strand"),
        )

    @classmethod
    def from_interval(cls, interval: IntervalCoordinates) -> GenomicInterval:
        """Create a value object from any protocol-compatible interval."""
        return cls.from_raw(
            chromosome=interval.chromosome,
            start=interval.minpos,
            end=interval.maxpos,
            strand=interval.strand,
        )

    @property
    def minpos(self) -> int:
        return self.start

    @property
    def maxpos(self) -> int:
        return self.end

    def overlaps(self, other: IntervalCoordinates) -> bool:
        if self.chromosome != other.chromosome:
            return False
        return self.minpos <= other.maxpos and other.minpos <= self.maxpos

    def contains(self, other: IntervalCoordinates) -> bool:
        if self.chromosome != other.chromosome:
            return False
        return self.minpos <= other.minpos and self.maxpos >= other.maxpos

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1


__all__ = ["GenomicInterval", "IntervalCoordinates"]
