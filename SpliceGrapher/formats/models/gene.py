"""Gene entities and gene-oriented helpers for gene-model domains."""

from __future__ import annotations

from collections.abc import Callable, Iterable
from dataclasses import InitVar, dataclass, field
from itertools import chain
from operator import attrgetter
from typing import TypeAlias

from SpliceGrapher.core.enums import RecordType

from .constants import ISOFORM_TYPE, AttrMap, normalize_attributes
from .features import BaseFeature, Exon, TranscriptRegion
from .locus import Locus
from .transcript import Transcript

GeneFilter: TypeAlias = Callable[["Gene"], bool]


def gene_type_filter(g: Gene) -> bool:
    """Convenience filter for getting only 'gene' records."""

    return g.feature_type == RecordType.GENE


def default_gene_filter(g: Gene) -> bool:
    """Default function for filtering genes from a list."""

    return True


@dataclass(slots=True, eq=True)
class Gene:
    id: str
    note: str | None
    start_pos: InitVar[int]
    end_pos: InitVar[int]
    chromosome: str
    strand: str
    name: str | None = None
    attr: InitVar[AttrMap | None] = None
    is_pseudo: bool = field(default=False, compare=False)
    minpos: int = field(init=False)
    maxpos: int = field(init=False)
    attributes: dict[str, str] = field(init=False, default_factory=dict, compare=False)
    transcripts: dict[str, Transcript] = field(default_factory=dict, compare=False)
    exons: list[Exon] = field(default_factory=list, compare=False)
    cds: list[TranscriptRegion] = field(default_factory=list, compare=False)
    exon_map: dict[tuple[int, int], Exon] = field(default_factory=dict, compare=False)
    cds_map: dict[tuple[str | RecordType, int, int], TranscriptRegion] = field(
        default_factory=dict,
        compare=False,
    )
    start_codon_map: dict[str, tuple[int, int] | None] = field(default_factory=dict, compare=False)
    end_codon_map: dict[str, tuple[int, int] | None] = field(default_factory=dict, compare=False)
    features: list[BaseFeature] = field(default_factory=list, compare=False)

    def __post_init__(
        self,
        start_pos: int,
        end_pos: int,
        attr: AttrMap | None,
    ) -> None:
        if self.strand not in {"+", "-", "."}:
            raise ValueError(f"Unknown strand '{self.strand}'")
        self.chromosome = self.chromosome.lower()
        self.minpos = min(start_pos, end_pos)
        self.maxpos = max(start_pos, end_pos)
        self.attributes = normalize_attributes(attr)
        if self.name is None:
            self.name = self.id

    @property
    def feature_type(self) -> RecordType:
        return RecordType.PSEUDOGENE if self.is_pseudo else RecordType.GENE

    @property
    def locus(self) -> Locus:
        return Locus.create(self.chromosome, self.minpos, self.maxpos, self.strand)

    def __len__(self) -> int:
        return len(self.locus)

    def contains(self, pos: int, strand: str) -> bool:
        return self.locus.contains(pos, strand)

    def start(self) -> int:
        return self.locus.start()

    def end(self) -> int:
        return self.locus.end()

    def _iter_transcripts(self) -> Iterable[Transcript]:
        yield from self.transcripts.values()

    def _validate_transcript_membership(self, transcript: Transcript) -> None:
        if transcript.strand != self.strand:
            raise ValueError(
                f"Transcript strand '{transcript.strand}' does not match gene strand "
                f"'{self.strand}' for {self.id}"
            )
        if transcript.chromosome != self.chromosome:
            raise ValueError(
                f"Transcript chromosome '{transcript.chromosome}' does not match gene "
                f"chromosome '{self.chromosome}' for {self.id}"
            )

    def _validate_region_membership(self, region: TranscriptRegion) -> None:
        if region.strand != self.strand:
            raise ValueError(
                f"Region strand '{region.strand}' does not match gene strand "
                f"'{self.strand}' for {self.id}"
            )
        if region.chromosome != self.chromosome:
            raise ValueError(
                f"Region chromosome '{region.chromosome}' does not match gene chromosome "
                f"'{self.chromosome}' for {self.id}"
            )

    def _register_gene_exon(self, new_exon: Exon) -> tuple[Exon, bool]:
        self._validate_region_membership(new_exon)
        pos_tuple = (new_exon.minpos, new_exon.maxpos)
        existing = self.exon_map.get(pos_tuple)
        if existing is None:
            self.exons.append(new_exon)
            self.exon_map[pos_tuple] = new_exon
            self.minpos = min(self.minpos, new_exon.minpos)
            self.maxpos = max(self.maxpos, new_exon.maxpos)
            return new_exon, True
        existing.attributes = new_exon.attributes | existing.attributes
        return existing, False

    def _register_gene_cds(self, new_cds: TranscriptRegion) -> tuple[TranscriptRegion, bool]:
        self._validate_region_membership(new_cds)
        cds_tuple = (new_cds.feature_type, new_cds.minpos, new_cds.maxpos)
        existing = self.cds_map.get(cds_tuple)
        if existing is None:
            self.cds.append(new_cds)
            self.cds_map[cds_tuple] = new_cds
            self.minpos = min(self.minpos, new_cds.minpos)
            self.maxpos = max(self.maxpos, new_cds.maxpos)
            return new_cds, True
        existing.attributes = new_cds.attributes | existing.attributes
        return existing, False

    def add_transcript(self, transcript: Transcript) -> Transcript:
        self._validate_transcript_membership(transcript)
        transcript.set_parent(self.id)
        existing = self.transcripts.get(transcript.id)
        if existing is None:
            self.transcripts[transcript.id] = transcript
            self.minpos = min(self.minpos, transcript.minpos)
            self.maxpos = max(self.maxpos, transcript.maxpos)
            for exon in transcript.exons:
                self._register_gene_exon(exon)
            for cds in transcript.cds:
                self._register_gene_cds(cds)
            self.start_codon_map[transcript.id] = transcript.start_codon()
            self.end_codon_map[transcript.id] = transcript.end_codon()
            return transcript

        for exon in transcript.exons:
            gene_exon, _ = self._register_gene_exon(exon)
            existing.add_exon(gene_exon)
        for cds in transcript.cds:
            gene_cds, _ = self._register_gene_cds(cds)
            existing.add_cds(gene_cds)
        for feature in transcript.features:
            existing.add_feature(feature)
        existing.minpos = min(existing.minpos, transcript.minpos)
        existing.maxpos = max(existing.maxpos, transcript.maxpos)
        existing.attributes.update(transcript.attributes)
        if existing.feature_type == ISOFORM_TYPE and transcript.feature_type != ISOFORM_TYPE:
            existing.feature_type = transcript.feature_type
        self.minpos = min(self.minpos, existing.minpos)
        self.maxpos = max(self.maxpos, existing.maxpos)
        self.start_codon_map[existing.id] = existing.start_codon()
        self.end_codon_map[existing.id] = existing.end_codon()
        return existing

    def add_exon(self, transcript: Transcript, new_exon: Exon) -> bool:
        stored_transcript = self.add_transcript(transcript)
        exon, result = self._register_gene_exon(new_exon)
        stored_transcript.add_exon(exon)
        return result

    def add_cds(self, transcript: Transcript, new_cds: TranscriptRegion) -> bool:
        stored_transcript = self.add_transcript(transcript)
        cds, result = self._register_gene_cds(new_cds)
        stored_transcript.add_cds(cds)
        self.start_codon_map[stored_transcript.id] = stored_transcript.start_codon()
        self.end_codon_map[stored_transcript.id] = stored_transcript.end_codon()
        return result

    def add_feature(self, feature: BaseFeature) -> None:
        if feature.strand != self.strand:
            raise ValueError(
                f"ERROR: feature strand '{feature.strand}' does not match gene strand "
                f"'{self.strand}'"
            )
        if feature.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: feature chromosome '{feature.chromosome}' does not match gene "
                f"chromosome '{self.chromosome}'"
            )

        feature.set_parent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def _collect_transcript_sites(
        self,
        site_fn: Callable[[Transcript], list[int]],
    ) -> list[int]:
        sites = chain.from_iterable(site_fn(transcript) for transcript in self._iter_transcripts())
        return sorted(set(sites), reverse=(self.strand == "-"))

    def acceptor_list(self) -> list[int]:
        return self._collect_transcript_sites(lambda t: t.acceptor_list())

    def donor_list(self) -> list[int]:
        return self._collect_transcript_sites(lambda t: t.donor_list())

    def end_codons(self) -> dict[str, tuple[int, int] | None]:
        return self.end_codon_map

    def get_feature_list(self, feature_type: str | RecordType) -> list[BaseFeature]:
        return [feature for feature in self.features if feature.feature_type == feature_type]

    def get_introns(self) -> set[tuple[int, int]]:
        introns: set[tuple[int, int]] = set()
        for transcript in self._iter_transcripts():
            exons = transcript.sorted_exons()
            for i in range(1, len(exons)):
                introns.add((exons[i - 1].end(), exons[i].start()))
        return introns

    def get_transcript(self, transcript_id: str) -> Transcript:
        return self.transcripts[transcript_id]

    def get_junctions(self) -> set[tuple[int, int]]:
        junctions: set[tuple[int, int]] = set()
        for transcript in self._iter_transcripts():
            exons = transcript.sorted_exons()
            for i in range(1, len(exons)):
                junctions.add((exons[i - 1].donor(), exons[i].acceptor()))
        return junctions

    def is_single_exon(self) -> bool:
        if not self.transcripts:
            return False
        return all(len(transcript.sorted_exons()) == 1 for transcript in self._iter_transcripts())

    def sorted_exons(self) -> list[Exon]:
        all_exons = chain.from_iterable(
            transcript.sorted_exons() for transcript in self._iter_transcripts()
        )
        by_coords: dict[tuple[int, int, str, str], Exon] = {}
        for exon in all_exons:
            key = (exon.minpos, exon.maxpos, exon.chromosome, exon.strand)
            if key not in by_coords:
                by_coords[key] = exon
        unique_exons = list(by_coords.values())
        return sorted(
            unique_exons,
            key=attrgetter("minpos", "maxpos"),
            reverse=(self.strand == "-"),
        )

    def start_codons(self) -> dict[str, tuple[int, int] | None]:
        return self.start_codon_map

    def detail_string(self) -> str:
        result = (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons) + len(self.cds)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}\n"
        )
        result += "\n  ".join(
            [transcript.detail_string() for transcript in self.transcripts.values()]
        )
        return result

    def __str__(self) -> str:
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.cds) + len(self.exons)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}"
        )


@dataclass(slots=True, eq=True)
class PseudoGene(Gene):
    is_pseudo: bool = field(default=True, init=False, compare=False)

    def detail_string(self) -> str:
        return self.__str__()
