from __future__ import annotations

from dataclasses import dataclass

import pytest

from SpliceGrapher.core.enums import Strand
from SpliceGrapher.core.types import GenomicInterval, IntervalCoordinates


@dataclass(slots=True)
class _LegacyInterval:
    chromosome: str
    minpos: int
    maxpos: int
    strand: str


def test_genomic_interval_orders_and_compares_as_value_object() -> None:
    left = GenomicInterval.from_raw("chr1", 5, 10, "+")
    right = GenomicInterval.from_raw("chr1", 11, 20, "+")
    same_as_left = GenomicInterval.from_raw("chr1", 5, 10, "+")

    assert left < right
    assert left == same_as_left
    assert len({left, same_as_left, right}) == 2


def test_genomic_interval_rejects_invalid_coordinates() -> None:
    with pytest.raises(ValueError, match=">= 1"):
        GenomicInterval.from_raw("chr1", 0, 10, "+")
    with pytest.raises(ValueError, match="start must be <= end"):
        GenomicInterval.from_raw("chr1", 11, 10, "+")


def test_genomic_interval_rejects_invalid_strand() -> None:
    with pytest.raises(ValueError, match="strand"):
        GenomicInterval.from_raw("chr1", 1, 10, "?")


def test_genomic_interval_interval_operations() -> None:
    parent = GenomicInterval.from_raw("chr1", 1, 20, Strand.PLUS)
    child = GenomicInterval.from_raw("chr1", 5, 10, Strand.PLUS)
    disjoint = GenomicInterval.from_raw("chr1", 21, 30, Strand.PLUS)

    assert parent.contains(child)
    assert parent.overlaps(child)
    assert not parent.overlaps(disjoint)


def test_genomic_interval_accepts_protocol_compatible_inputs() -> None:
    legacy = _LegacyInterval(chromosome="chr2", minpos=3, maxpos=8, strand="+")

    assert isinstance(legacy, IntervalCoordinates)

    interval = GenomicInterval.from_interval(legacy)
    assert interval.chromosome == "chr2"
    assert interval.minpos == 3
    assert interval.maxpos == 8
    assert interval.strand is Strand.PLUS
