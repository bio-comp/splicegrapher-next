"""Transcript entity for gene-model domains."""

from __future__ import annotations

from dataclasses import InitVar, dataclass, field, replace
from itertools import pairwise
from operator import attrgetter

from SpliceGrapher.core.enums import RecordType

from .constants import AttrMap, normalize_attributes
from .features import BaseFeature, Exon, TranscriptRegion
from .index import _SpliceSiteLike, feature_sort_key
from .locus import Locus


@dataclass(slots=True, eq=True)
class Transcript:
    id: str
    start_pos: InitVar[int]
    end_pos: InitVar[int]
    chromosome_value: InitVar[str]
    strand_value: InitVar[str]
    attr: InitVar[AttrMap | None] = None
    feature_type: str | RecordType = field(default=RecordType.MRNA, compare=False)
    parent: str | None = field(default=None, compare=False)
    locus: Locus = field(init=False)
    attributes: dict[str, str] = field(init=False, default_factory=dict, compare=False)
    features: list[BaseFeature] = field(default_factory=list, compare=False)
    exons: list[Exon] = field(default_factory=list, compare=False)
    exon_map: dict[tuple[int, int], Exon] = field(default_factory=dict, compare=False)
    cds: list[TranscriptRegion] = field(default_factory=list, compare=False)
    cds_map: dict[tuple[str | RecordType, int, int], TranscriptRegion] = field(
        default_factory=dict,
        compare=False,
    )
    start_codon_pos: tuple[int, int] | None = field(default=None, compare=False)
    end_codon_pos: tuple[int, int] | None = field(default=None, compare=False)

    def __post_init__(
        self,
        start_pos: int,
        end_pos: int,
        chromosome_value: str,
        strand_value: str,
        attr: AttrMap | None,
    ) -> None:
        self.locus = Locus.create(chromosome_value, start_pos, end_pos, strand_value)
        self.attributes = normalize_attributes(attr)

    @property
    def chromosome(self) -> str:
        return self.locus.chromosome

    @property
    def strand(self) -> str:
        return self.locus.strand

    @property
    def minpos(self) -> int:
        return self.locus.minpos

    @minpos.setter
    def minpos(self, value: int) -> None:
        lo = min(value, self.locus.maxpos)
        hi = max(value, self.locus.maxpos)
        self.locus = replace(self.locus, minpos=lo, maxpos=hi)

    @property
    def maxpos(self) -> int:
        return self.locus.maxpos

    @maxpos.setter
    def maxpos(self, value: int) -> None:
        lo = min(self.locus.minpos, value)
        hi = max(self.locus.minpos, value)
        self.locus = replace(self.locus, minpos=lo, maxpos=hi)

    @property
    def record_type(self) -> str | RecordType:
        return self.feature_type

    def set_parent(self, parent_id: str) -> None:
        self.parent = parent_id

    def start(self) -> int:
        return self.locus.start()

    def end(self) -> int:
        return self.locus.end()

    def contains(self, pos: int, strand: str) -> bool:
        return self.locus.contains(pos, strand)

    def __len__(self) -> int:
        return len(self.locus)

    def add_exon(self, exon: Exon) -> bool:
        if exon.strand != self.strand:
            raise ValueError(
                f"Exon strand '{exon.strand}' does not match isoform strand '{self.strand}' "
                f"for {self.id}"
            )
        if exon.chromosome != self.chromosome:
            raise ValueError(
                f"Exon chromosome '{exon.chromosome}' does not match isoform chromosome "
                f"'{self.chromosome}' for {self.id}"
            )

        exon_tuple = (exon.minpos, exon.maxpos)
        if exon_tuple in self.exon_map:
            existing = self.exon_map[exon_tuple]
            existing.attributes = exon.attributes | existing.attributes
            return False

        self.exons.append(exon)
        self.exon_map[exon_tuple] = exon
        self.minpos = min(self.minpos, exon.minpos)
        self.maxpos = max(self.maxpos, exon.maxpos)
        return True

    def add_cds(self, cds: TranscriptRegion) -> bool:
        if cds.strand != self.strand:
            raise ValueError(
                f"ERROR: CDS strand '{cds.strand}' does not match gene strand '{self.strand}'"
            )
        if cds.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: CDS chromosome '{cds.chromosome}' does not match "
                f"gene chromosome '{self.chromosome}'"
            )

        cds_tuple = (cds.feature_type, cds.minpos, cds.maxpos)
        if cds_tuple in self.cds_map:
            existing = self.cds_map[cds_tuple]
            existing.attributes = cds.attributes | existing.attributes
            return False
        self.cds_map[cds_tuple] = cds
        self.cds.append(cds)
        self.minpos = min(self.minpos, cds.minpos)
        self.maxpos = max(self.maxpos, cds.maxpos)
        return True

    def add_feature(self, feature: BaseFeature) -> None:
        if feature.strand != self.strand:
            raise ValueError(
                f"ERROR: feature strand '{feature.strand}' does not match form strand "
                f"'{self.strand}'"
            )
        if feature.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: feature chromosome '{feature.chromosome}' does not match form "
                f"chromosome '{self.chromosome}'"
            )

        feature.set_parent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def _ordered_splice_regions(self) -> list[_SpliceSiteLike]:
        regions: list[_SpliceSiteLike] = list(self.cds) if self.cds else list(self.exons)
        return sorted(
            regions,
            key=attrgetter("minpos", "maxpos"),
            reverse=(self.strand == "-"),
        )

    def _ordered_cds_regions(self) -> list[TranscriptRegion]:
        return sorted(
            self.cds,
            key=attrgetter("minpos", "maxpos"),
            reverse=(self.strand == "-"),
        )

    def _cds_start_bounds(self, region: TranscriptRegion) -> tuple[int, int]:
        return (
            (region.minpos, region.minpos + 2)
            if self.strand == "+"
            else (region.maxpos - 2, region.maxpos)
        )

    def _cds_end_bounds(self, region: TranscriptRegion) -> tuple[int, int]:
        return (
            (region.maxpos - 2, region.maxpos)
            if self.strand == "+"
            else (region.minpos, region.minpos + 2)
        )

    def acceptor_list(self) -> list[int]:
        return [region.acceptor() for region in self._ordered_splice_regions()[1:]]

    def donor_list(self) -> list[int]:
        return [region.donor() for region in self._ordered_splice_regions()[:-1]]

    def end_codon(self) -> tuple[int, int] | None:
        if not self.end_codon_pos:
            self.find_codons()
        return self.end_codon_pos

    def find_codons(self) -> None:
        if self.end_codon_pos and self.start_codon_pos:
            return
        ordered_cds = self._ordered_cds_regions()
        if not ordered_cds:
            return
        for prev, c in pairwise(ordered_cds):
            if (
                not self.start_codon_pos
                and prev.feature_type == RecordType.FIVE_PRIME_UTR
                and c.feature_type == RecordType.CDS
            ):
                self.start_codon_pos = self._cds_start_bounds(c)
            elif (
                not self.end_codon_pos
                and prev.feature_type == RecordType.CDS
                and c.feature_type == RecordType.THREE_PRIME_UTR
            ):
                self.end_codon_pos = self._cds_end_bounds(prev)
            if self.start_codon_pos and self.end_codon_pos:
                break

    def infer_codons(self) -> None:
        ordered_cds = self._ordered_cds_regions()
        if not ordered_cds:
            return

        self.find_codons()
        if not self.start_codon_pos:
            self.start_codon_pos = self._cds_start_bounds(ordered_cds[0])

        if not self.end_codon_pos:
            self.end_codon_pos = self._cds_end_bounds(ordered_cds[-1])

    def get_utrs(self) -> list[TranscriptRegion]:
        return [
            c
            for c in self.cds
            if c.feature_type in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]
        ]

    def get_feature_list(self, feature_type: str | RecordType) -> list[BaseFeature]:
        return [f for f in self.features if f.feature_type == feature_type]

    def sorted_cds(self) -> list[TranscriptRegion]:
        return self._ordered_cds_regions()

    def sorted_exons(self, *, minintron: int = 2) -> list[Exon]:
        if self.exons:
            return sorted(self.exons, key=feature_sort_key, reverse=(self.strand == "-"))

        cds_list = sorted(self.cds, key=attrgetter("minpos", "maxpos"))
        result: list[Exon] = []
        for cds in cds_list:
            gap = cds.minpos - result[-1].maxpos if result else None
            if gap is not None and 0 <= gap < minintron:
                prev = result[-1]
                result[-1] = Exon(
                    prev.minpos,
                    cds.maxpos,
                    self.chromosome,
                    self.strand,
                    attr=prev.attributes,
                )
            else:
                result.append(Exon(cds.minpos, cds.maxpos, self.chromosome, self.strand))
        if self.strand == "-":
            result.reverse()
        for exon in result:
            exon.feature_type = RecordType.CDS
        return result

    def sorted_introns(self) -> list[tuple[int, int]]:
        result: list[tuple[int, int]] = []
        exons = self.sorted_exons()
        if len(exons) < 2:
            return result
        prev = exons[0]
        for exon in exons[1:]:
            result.append((prev.donor(), exon.acceptor()))
            prev = exon
        return result

    def start_codon(self) -> tuple[int, int] | None:
        if not self.start_codon_pos:
            self.find_codons()
        return self.start_codon_pos

    def exon_string(self) -> str:
        return ",".join([str(exon) for exon in self.exons])

    def detail_string(self) -> str:
        return (
            f"Isoform {self.id}\nStart: {self.start()}; End: {self.end()}; "
            f"Strand: {self.strand}\nExons: [{self.exon_string()}]\n"
        )

    def __str__(self) -> str:
        feature_count = len(self.exons) + len(self.cds)
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {feature_count} exons/cds, "
            f"range {self.minpos} to {self.maxpos}"
        )
