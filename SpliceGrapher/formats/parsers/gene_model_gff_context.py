"""Parser state and protocols for ``GeneModel`` GFF loading."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Protocol

import structlog

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.core.enums import RecordType

if TYPE_CHECKING:
    from SpliceGrapher.formats.parsers.gene_model_gff_records import ParsedRecord

LOGGER = structlog.get_logger(__name__)
PARENT_CACHE_MAX_ENTRIES = 200_000
PARENT_CANDIDATE_MAX_ENTRIES = 200_000


class GeneModelLike(Protocol):
    model: dict[str, dict[str, model_domain.Gene]]
    mrna_forms: dict[str, dict[str, model_domain.Transcript]]
    mrna_gene: dict[str, dict[str, model_domain.Gene]]
    all_genes: dict[str, model_domain.Gene]
    found_types: dict[RecordType, bool]
    all_chr: dict[str, model_domain.Chromosome]
    chromosome_index: dict[str, model_domain.ChromosomeGeneIndex]

    def add_chromosome(self, start: int, end: int, name: str) -> None: ...

    def add_gene(self, gene: model_domain.Gene) -> None: ...

    def get_parent(
        self,
        s: str,
        chrom: str,
        search_genes: bool = True,
        search_mrna: bool = True,
    ) -> model_domain.Gene | model_domain.Transcript | None: ...

    def get_mrna_parent(
        self,
        chrom: str,
        mrna_id: str,
    ) -> model_domain.Gene | None: ...


@dataclass(slots=True)
class ParseStats:
    gene_count: int = 0
    exon_count: int = 0
    iso_count: int = 0
    mrna_count: int = 0
    cds_count: int = 0
    orphan_mrna_count: int = 0
    bad_lines: int = 0
    parent_cache_clear_count: int = 0
    parent_candidate_clear_count: int = 0


@dataclass(slots=True)
class ParseContext:
    model: GeneModelLike
    require_notes: bool
    verbose: bool
    ignore_errors: bool
    chromosomes: set[str] | None = None
    parent_cache_max_entries: int = PARENT_CACHE_MAX_ENTRIES
    parent_candidate_max_entries: int = PARENT_CANDIDATE_MAX_ENTRIES
    stats: ParseStats = field(default_factory=ParseStats)
    gene_alias: dict[str, str] = field(default_factory=dict)
    parent_cache: dict[
        tuple[str, str, bool, bool],
        model_domain.Gene | model_domain.Transcript,
    ] = field(default_factory=dict)
    parent_candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)
    pending_exons: dict[tuple[str, str], list[ParsedRecord]] = field(default_factory=dict)
    pending_regions: dict[tuple[str, str], list[ParsedRecord]] = field(default_factory=dict)

    def fail(self, message: str) -> None:
        if self.ignore_errors:
            return
        raise RuntimeError(message)

    def write_verbose(self, message: str) -> None:
        if self.verbose:
            LOGGER.info("gene_model_gff_verbose", message=message)

    def report_bad_line(self, line_no: int) -> None:
        self.stats.bad_lines += 1
        self.write_verbose(
            f"line {line_no}: invalid GFF format (not enough columns); file may be corrupt"
        )
        if self.stats.bad_lines >= model_domain.MAX_BAD_LINES:
            LOGGER.error("gene_model_gff_invalid_input", line=line_no)
            raise ValueError("Invalid GFF input file")

    def cache_parent_candidate(self, parent_string: str, candidates: tuple[str, ...]) -> None:
        if (
            self.parent_candidate_max_entries > 0
            and parent_string not in self.parent_candidates
            and len(self.parent_candidates) >= self.parent_candidate_max_entries
        ):
            self.parent_candidates.clear()
            self.stats.parent_candidate_clear_count += 1
        self.parent_candidates[parent_string] = candidates

    def cache_parent(
        self,
        cache_key: tuple[str, str, bool, bool],
        parent_record: model_domain.Gene | model_domain.Transcript,
    ) -> None:
        if (
            self.parent_cache_max_entries > 0
            and cache_key not in self.parent_cache
            and len(self.parent_cache) >= self.parent_cache_max_entries
        ):
            self.parent_cache.clear()
            self.stats.parent_cache_clear_count += 1
        self.parent_cache[cache_key] = parent_record


__all__ = [
    "GeneModelLike",
    "LOGGER",
    "PARENT_CACHE_MAX_ENTRIES",
    "PARENT_CANDIDATE_MAX_ENTRIES",
    "ParseContext",
    "ParseStats",
]
