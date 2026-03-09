"""Gene-model interval/index helpers and sort utilities.

Split from `SpliceGrapher.formats.models` to keep the public import surface
stable while shrinking the domain monolith.
"""

from __future__ import annotations

from collections.abc import Iterable
from dataclasses import dataclass
from typing import TYPE_CHECKING, Protocol

from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)

if TYPE_CHECKING:
    from .domain import Gene


@dataclass(slots=True)
class IntervalQuery:
    """Simple interval wrapper for index overlap queries."""

    minpos: int
    maxpos: int


@dataclass(slots=True)
class ChromosomeGeneIndex:
    """Per-chromosome interval/search indexes keyed by strand."""

    by_strand: dict[str, list[Gene]]
    interval_by_strand: dict[str, InMemoryIntervalIndex[Gene]]

    @classmethod
    def build(cls, genes: Iterable[Gene]) -> ChromosomeGeneIndex:
        by_strand: dict[str, list[Gene]] = {"-": [], "+": [], ".": []}
        for gene in sorted(genes, key=gene_sort_key):
            by_strand.setdefault(gene.strand, []).append(gene)

        interval_by_strand = {
            strand: InMemoryIntervalIndex(strand_genes)
            for strand, strand_genes in by_strand.items()
            if strand_genes
        }
        return cls(by_strand=by_strand, interval_by_strand=interval_by_strand)

    def find_gene(self, start_pos: int, end_pos: int, strand: str) -> Gene | None:
        strand_index = self.interval_by_strand.get(strand)
        if strand_index is None:
            return None
        query = IntervalQuery(min(start_pos, end_pos), max(start_pos, end_pos))
        for gene in strand_index.overlaps(query):
            if gene.contains(start_pos, strand) or gene.contains(end_pos, strand):
                return gene
        return None

    def find_genes_overlapping(self, start_pos: int, end_pos: int, strand: str) -> list[Gene]:
        strand_index = self.interval_by_strand.get(strand)
        if strand_index is None:
            return []
        query = IntervalQuery(min(start_pos, end_pos), max(start_pos, end_pos))
        return list(strand_index.overlaps(query))

    def genes(self, strand: str | None = None) -> list[Gene]:
        if strand is None:
            return [gene for strand_genes in self.by_strand.values() for gene in strand_genes]
        return list(self.by_strand.get(strand, []))


@dataclass(slots=True)
class Chromosome:
    """Class that encapsulates a chromosome GFF record."""

    minpos: int
    maxpos: int
    name: str

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def __str__(self) -> str:
        return f"{self.name}: {self.minpos}-{self.maxpos}"

    def contains(self, pos: int) -> bool:
        return self.minpos <= pos <= self.maxpos

    def end(self) -> int:
        return self.maxpos

    def start(self) -> int:
        return self.minpos

    def update(self, feature: _IntervalLike) -> None:
        """
        Many species do not have 'chromosome' entries in their annotations,
        so we must infer chromosome boundaries from observed features.
        """
        if feature.chromosome.lower() != self.name.lower():
            raise ValueError(f"Cannot use feature from {feature.chromosome} to update {self.name}")
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)


def feature_cmp(a: _IntervalLike, b: _IntervalLike) -> int:
    """
    General comparison function for sorting any features that have 'minpos'
    and 'maxpos' attributes.
    """
    if a.minpos == b.minpos:
        return a.maxpos - b.maxpos
    else:
        return a.minpos - b.minpos


def feature_sort_key(feature: _IntervalLike) -> tuple[int, int]:
    """Sort features by genomic interval."""
    return (feature.minpos, feature.maxpos)


def gene_sort_key(gene: Gene) -> tuple[int, int, str]:
    """Sort genes by interval, then id for deterministic ties."""
    return (gene.minpos, gene.maxpos, gene.id)


def gtf_feature_sort_key(feature: _IntervalLike) -> tuple[int, int]:
    """Sort key for GTF emission using the genomic-ascending compatibility policy."""
    return feature_sort_key(feature)


def feature_overlaps(a: _IntervalLike | None, b: _IntervalLike | None) -> bool:
    """
    General function for determining whether feature 'a' and feature 'b' overlap.
    """
    if not (a and b):
        return False
    return intervals_overlap(a, b)


def feature_contains(a: _IntervalLike | None, b: _IntervalLike | None) -> bool:
    """
    General function for determining whether feature 'a' contains
    feature 'b'.  Note that both features must have 'minpos'
    and 'maxpos' attributes.
    """
    if not (a and b):
        return False
    return interval_contains(a, b)


class _IntervalLike(Protocol):
    minpos: int
    maxpos: int
    chromosome: str
    strand: str


class _SpliceSiteLike(_IntervalLike, Protocol):
    def donor(self) -> int: ...

    def acceptor(self) -> int: ...
