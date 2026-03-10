"""Canonical splice-junction model and parser for SGN depth records."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import cast

from SpliceGrapher.core.enums import JunctionCode, ShortReadCode

_JUNCTION_CODES = frozenset(code.value for code in JunctionCode)


@dataclass(eq=False)
class SpliceJunction:
    """Splice-junction record used by alignment and depth readers."""

    chromosome: str
    start: int
    end: int
    anchors: list[int] | tuple[int, int]
    sj_code: str
    strand: str
    count: int = 1
    code: ShortReadCode = field(init=False, default=ShortReadCode.JUNCTION)

    def __post_init__(self) -> None:
        self.chromosome = self.chromosome.lower()
        if self.start > self.end:
            self.start, self.end = self.end, self.start
        anchor_values = [int(value) for value in self.anchors]
        if len(anchor_values) != 2:
            raise ValueError("SpliceJunction anchors must contain exactly two values")
        self.anchors = anchor_values

    @property
    def p1(self) -> int:
        return self.start

    @property
    def p2(self) -> int:
        return self.end

    @property
    def minpos(self) -> int:
        return self.start

    @property
    def maxpos(self) -> int:
        return self.end

    @property
    def accval(self) -> int:
        return self.end if self.strand == "+" else self.start

    @property
    def donval(self) -> int:
        return self.start if self.strand == "+" else self.end

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SpliceJunction):
            return NotImplemented
        return (
            self.chromosome == other.chromosome
            and self.strand == other.strand
            and self.start == other.start
            and self.end == other.end
        )

    def __lt__(self, other: SpliceJunction) -> bool:
        if self.chromosome != other.chromosome:
            return self.chromosome < other.chromosome
        if self.start != other.start:
            return self.start < other.start
        return self.end < other.end

    def __hash__(self) -> int:
        return hash((self.chromosome, self.start, self.end, self.strand))

    def __repr__(self) -> str:
        return (
            "Junction: "
            f"{self.count} reads on {self.chromosome} between donor {self.donval} "
            f"and acceptor {self.accval} ({self.strand})"
        )

    def __str__(self) -> str:
        return f"{self.chromosome} don={self.donval}, acc={self.accval} ({self.strand})"

    def acceptor(self, strand: str | None = None) -> int:
        """Return the acceptor site for the junction based on strand."""
        if strand is None:
            return self.accval
        return self.end if strand == "+" else self.start

    def donor(self, strand: str | None = None) -> int:
        """Return the donor site for the junction based on strand."""
        if strand is None:
            return self.donval
        return self.start if strand == "+" else self.end

    def min_anchor(self) -> int:
        """Return the smaller of the two anchor values."""
        return min(self.anchors)

    def to_string(self) -> str:
        """Return tab-delimited SGN depth representation."""
        return (
            f"{ShortReadCode.JUNCTION.value}\t{self.chromosome}\t{self.strand}\t"
            f"{self.start}\t{self.end}\t{self.anchors[0]}\t{self.anchors[1]}"
            f"\t{self.sj_code}\t{self.count}"
        )

    def update(self, other: SpliceJunction) -> None:
        """Merge another equivalent junction into this record."""
        if self != other:
            raise ValueError(f"Splice sites do not match: {other} <--> {self}")
        self.count += other.count
        self_anchors = cast(list[int], self.anchors)
        other_anchors = cast(list[int], other.anchors)
        self_anchors[0] = max(self_anchors[0], other_anchors[0])
        self_anchors[1] = max(self_anchors[1], other_anchors[1])


def parse_junction_record(record: str) -> SpliceJunction:
    """Parse one legacy ``J`` depth record into a ``SpliceJunction``."""
    try:
        code, chrom, strand, p1, p2, a1, a2, sj_code, count = record.rstrip("\n").split("\t")
    except ValueError as exc:
        raise ValueError(f"Invalid SpliceJunction record format: {record!r}") from exc

    if code != ShortReadCode.JUNCTION.value:
        raise ValueError(f"Invalid SpliceJunction code: {code}")

    if sj_code not in _JUNCTION_CODES:
        raise ValueError(f"Invalid SpliceJunction type: {sj_code}")

    return SpliceJunction(
        chromosome=chrom,
        start=int(p1),
        end=int(p2),
        anchors=[int(a1), int(a2)],
        sj_code=sj_code,
        strand=strand,
        count=int(count),
    )


__all__ = ["SpliceJunction", "parse_junction_record"]
