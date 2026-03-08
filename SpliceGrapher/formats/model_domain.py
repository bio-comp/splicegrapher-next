"""Gene-model domain entities and shared constants.

Split from `SpliceGrapher.formats.models` to keep the public import surface
stable while shrinking the domain monolith.
"""

from __future__ import annotations

# TRYING A NEW WAY TO IDENTIFY GENES:
from collections.abc import Callable, Iterable, Mapping
from dataclasses import InitVar, dataclass, field, replace
from itertools import chain, pairwise
from operator import attrgetter
from typing import TypeVar

from SpliceGrapher.core.enums import AttrKey, RecordType, Strand
from SpliceGrapher.formats.model_index import _SpliceSiteLike, feature_sort_key

KNOWN_RECTYPES = [
    RecordType.CDS,
    RecordType.CHROMOSOME,
    RecordType.EXON,
    RecordType.FIVE_PRIME_UTR,
    RecordType.GENE,
    RecordType.INTRON,
    RecordType.MRNA,
    RecordType.MRNA_TE_GENE,
    RecordType.NONUNIQUE,
    RecordType.PROTEIN,
    RecordType.CDS_PREDICTED,
    RecordType.PREDICTED_GENE,
    RecordType.PSEUDOGENE,
    RecordType.PSEUDOGENIC_EXON,
    RecordType.PSEUDOGENIC_TRANSCRIPT,
    RecordType.THREE_PRIME_UTR,
    RecordType.TRANS_ELE_GENE,
]
IGNORE_RECTYPES = {
    RecordType.PROTEIN,
    RecordType.INTRON,
    RecordType.MRNA_TE_GENE,
    RecordType.TRANS_ELE_GENE,
    RecordType.NONUNIQUE,
}
CDS_TYPES = {RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR, RecordType.CDS}

# Special GTF types for conversions
GTF_GENE_ID = "gene_id"
GTF_GENE_NAME = "gene_name"
GTF_TRANSCRIPT = "transcript_id"
GTF_TRANSNAME = "transcript_name"
GTF_SOURCE = "gene_biotype"
GTF_EXON_ID = "exon_number"
GTF_PROTEIN_ID = "protein_id"

# Record type map allows mapping unusual names to known types:
RECTYPE_MAP = {k: k for k in KNOWN_RECTYPES}
RECTYPE_MAP[RecordType.PREDICTED_GENE] = RecordType.GENE
RECTYPE_MAP[RecordType.CDS_PREDICTED] = RecordType.CDS

# Virtual types (don't appear in GFF):
ISOFORM_TYPE = "isoform"

# Annotation fields:
ID_FIELD = AttrKey.ID
NAME_FIELD = AttrKey.NAME
NOTE_FIELD = AttrKey.NOTE
PARENT_FIELD = AttrKey.PARENT

# There seems to be no consensus about how these tags are used in UCSC, ENSEMBL
# and other forms of gene models, so we must try each kind:
POSSIBLE_GENE_FIELDS = [PARENT_FIELD, GTF_GENE_ID, GTF_GENE_NAME]
POSSIBLE_FORM_FIELDS = [PARENT_FIELD, GTF_TRANSCRIPT]

# ID for GFF output:
GFF_ID = "SpliceGrapher"
MAX_BAD_LINES = 3

# Some files contain '.' for an unassigned strand
VALID_STRANDS = {strand.value for strand in Strand}

FORM_DELIMITERS = [".", "-", "_", ","]

# GeneModel splice-site helpers assume canonical GT/AG-style 2nt dimer boundaries.
# This remains a compatibility contract for existing callers and fixtures.
SPLICE_DIMER_OFFSET = 2

# GTF serialization contract: emit records in genomic coordinate order for both
# strands (ascending minpos, then maxpos) to preserve historical output shape.
GTF_ORDER_POLICY = "genomic_ascending"

GeneFilter = Callable[["Gene"], bool]
GffRecordSource = str | list[str] | set[str] | tuple[str, ...]
AttrValue = str
AttrMap = Mapping[str | AttrKey, AttrValue]
ExtraAttrMap = Mapping[str, AttrValue] | Mapping[AttrKey, AttrValue]
AttrT = TypeVar("AttrT")


def gene_type_filter(g: Gene) -> bool:
    """Convenience filter for getting only 'gene' records."""
    return g.feature_type == RecordType.GENE


def default_gene_filter(g: Gene) -> bool:
    """Default function for filtering genes from a list."""
    return True


def cds_factory(
    rec_type: RecordType,
    start_pos: int,
    end_pos: int,
    chr_name: str,
    strand: str,
    attr: AttrMap | None = None,
) -> TranscriptRegion:
    """Simple factory method for creating CDS-type records."""
    attr = {} if attr is None else {str(key): str(value) for key, value in attr.items()}
    if rec_type == RecordType.CDS:
        return CDS(start_pos, end_pos, chr_name, strand, attr)
    elif rec_type == RecordType.FIVE_PRIME_UTR:
        return FpUtr(start_pos, end_pos, chr_name, strand, attr)
    elif rec_type == RecordType.THREE_PRIME_UTR:
        return TpUtr(start_pos, end_pos, chr_name, strand, attr)
    else:
        raise ValueError(f"Illegal CDS record type: {rec_type}")


def _normalize_attributes(attr: AttrMap | None) -> dict[str, str]:
    return {} if attr is None else {str(key): str(value) for key, value in attr.items()}


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
        self.attributes = _normalize_attributes(attr)

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
        default=RecordType.FIVE_PRIME_UTR, init=False, compare=False
    )


@dataclass(slots=True, eq=True)
class TpUtr(TranscriptRegion):
    feature_type: str | RecordType = field(
        default=RecordType.THREE_PRIME_UTR, init=False, compare=False
    )


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
        default_factory=dict, compare=False
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
        self.attributes = _normalize_attributes(attr)

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
        if self.strand not in VALID_STRANDS:
            raise ValueError(f"Unknown strand '{self.strand}'")
        self.chromosome = self.chromosome.lower()
        self.minpos = min(start_pos, end_pos)
        self.maxpos = max(start_pos, end_pos)
        self.attributes = _normalize_attributes(attr)
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
