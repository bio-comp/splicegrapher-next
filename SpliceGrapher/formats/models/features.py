"""Base feature and transcript-region entities for gene models."""

from __future__ import annotations

from dataclasses import InitVar, dataclass, field
from typing import TYPE_CHECKING

from SpliceGrapher.core.enums import RecordType

from .constants import VALID_STRANDS, AttrMap, normalize_attributes
from .locus import Locus

if TYPE_CHECKING:
    from .transcript import Transcript


@dataclass(slots=True, eq=True)
class BaseFeature:
    feature_type: str | RecordType = field(compare=False)
    start_pos: InitVar[int]
    end_pos: InitVar[int]
    chromosome: str
    strand: str
    attr: InitVar[AttrMap | None] = None
    parent: str | None = field(default=None, compare=False)
    minpos: int = field(init=False)
    maxpos: int = field(init=False)
    attributes: dict[str, str] = field(init=False, default_factory=dict, compare=False)

    def __post_init__(
        self,
        start_pos: int,
        end_pos: int,
        attr: AttrMap | None,
    ) -> None:
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"Unknown strand '{self.strand}'")
        self.chromosome = self.chromosome.lower()
        self.minpos = min(start_pos, end_pos)
        self.maxpos = max(start_pos, end_pos)
        self.attributes = normalize_attributes(attr)

    @property
    def locus(self) -> Locus:
        return Locus.create(self.chromosome, self.minpos, self.maxpos, self.strand)

    def start(self) -> int:
        return self.locus.start()

    def end(self) -> int:
        return self.locus.end()

    def contains(self, pos: int, strand: str) -> bool:
        return self.locus.contains(pos, strand)

    def donor(self) -> int:
        return self.locus.donor()

    def acceptor(self) -> int:
        return self.locus.acceptor()

    def detail_string(self) -> str:
        return (
            f"Feature: {self.feature_type}\nChromosome: {self.chromosome}\n"
            f"Start: {self.start()}; End: {self.end()}; Strand: '{self.strand}'"
        )

    def set_parent(self, parent_id: str) -> None:
        self.parent = parent_id

    def __len__(self) -> int:
        return len(self.locus)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
        )

    def __str__(self) -> str:
        return f"{self.chromosome} {self.feature_type}: {self.start()}-{self.end()} ({self.strand})"


@dataclass(slots=True, eq=True)
class TranscriptRegion(BaseFeature):
    def add_parent(self, isoform: Transcript) -> None:
        _ = isoform

    def __str__(self) -> str:
        return f"{self.feature_type} {self.start()}-{self.end()}({self.strand})"


@dataclass(slots=True, eq=True)
class Exon(TranscriptRegion):
    feature_type: str | RecordType = field(default=RecordType.EXON, init=False, compare=False)


@dataclass(slots=True, eq=True)
class CDS(TranscriptRegion):
    feature_type: str | RecordType = field(default=RecordType.CDS, init=False, compare=False)

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos, self.feature_type) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
            other.feature_type,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (
            self.feature_type == other.feature_type
            and self.minpos == other.minpos
            and self.maxpos == other.maxpos
            and self.strand == other.strand
            and self.chromosome == other.chromosome
        )

    def __str__(self) -> str:
        return f"{self.feature_type} {self.start()}-{self.end()}({self.strand})"


@dataclass(slots=True, eq=True)
class FpUtr(TranscriptRegion):
    feature_type: str | RecordType = field(
        default=RecordType.FIVE_PRIME_UTR,
        init=False,
        compare=False,
    )


@dataclass(slots=True, eq=True)
class TpUtr(TranscriptRegion):
    feature_type: str | RecordType = field(
        default=RecordType.THREE_PRIME_UTR,
        init=False,
        compare=False,
    )


def cds_factory(
    rec_type: RecordType,
    start_pos: int,
    end_pos: int,
    chr_name: str,
    strand: str,
    attr: AttrMap | None = None,
) -> TranscriptRegion:
    """Simple factory method for creating CDS-type records."""

    attr = normalize_attributes(attr)
    if rec_type == RecordType.CDS:
        return CDS(start_pos, end_pos, chr_name, strand, attr)
    if rec_type == RecordType.FIVE_PRIME_UTR:
        return FpUtr(start_pos, end_pos, chr_name, strand, attr)
    if rec_type == RecordType.THREE_PRIME_UTR:
        return TpUtr(start_pos, end_pos, chr_name, strand, attr)
    raise ValueError(f"Illegal CDS record type: {rec_type}")
