"""
Stores gene annotation information from a GFF3 annotation
file and provides methods for searching on the data.
"""

from __future__ import annotations

# TRYING A NEW WAY TO IDENTIFY GENES:
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import InitVar, dataclass, field
from typing import TypeVar

from SpliceGrapher.core.enums import AttrKey, RecordType, Strand
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)

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

    def genes(self, strand: str | None = None) -> list[Gene]:
        if strand is None:
            return [gene for strand_genes in self.by_strand.values() for gene in strand_genes]
        return list(self.by_strand.get(strand, []))


def gene_type_filter(g: Gene) -> bool:
    """Convenience filter for getting only 'gene' records."""
    return g.feature_type == RecordType.GENE


def default_gene_filter(g: Gene) -> bool:
    """Default function for filtering genes from a list."""
    return True


def dict_to_gff(d: Mapping[str, str]) -> str:
    """Returns a string representation of a dictionary based on the GFF3 annotation format."""
    return ";".join([f"{k}={v}" for k, v in sorted(d.items()) if k != "parent"])


def dict_to_gtf(d: Mapping[str, str]) -> str:
    """Returns a string representation of a dictionary based on the GTF annotation format."""
    return "; ".join([f'{k} "{v}"' for k, v in sorted(d.items()) if k != "parent"])


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

    def gff_string(self) -> str:
        # Example: Chr1    TAIR9   chromosome  1   30427671    .   .   .   ID=Chr1;Name=Chr1
        name_str = self.name.capitalize()
        return (
            f"{name_str}\t{GFF_ID}\tchromosome\t{self.start()}\t{self.end()}\t.\t.\t."
            f"\tID={name_str};Name={name_str}"
        )

    def start(self) -> int:
        return self.minpos

    def update(self, feature: BaseFeature) -> None:
        """
        Many species do not have 'chromosome' entries in their annotations,
        so we must infer chromosome boundaries from observed features.
        """
        if feature.chromosome.lower() != self.name.lower():
            raise ValueError(f"Cannot use feature from {feature.chromosome} to update {self.name}")
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)


def feature_cmp(a: BaseFeature, b: BaseFeature) -> int:
    """
    General comparison function for sorting any features that have 'minpos'
    and 'maxpos' attributes.
    """
    if a.minpos == b.minpos:
        return a.maxpos - b.maxpos
    else:
        return a.minpos - b.minpos


def feature_sort_key(feature: BaseFeature) -> tuple[int, int]:
    """Sort features by genomic interval."""
    return (feature.minpos, feature.maxpos)


def gene_sort_key(gene: Gene) -> tuple[int, int, str]:
    """Sort genes by interval, then id for deterministic ties."""
    return (gene.minpos, gene.maxpos, gene.id)


def gtf_feature_sort_key(feature: BaseFeature) -> tuple[int, int]:
    """Sort key for GTF emission using the genomic-ascending compatibility policy."""
    return feature_sort_key(feature)


def feature_overlaps(a: BaseFeature | None, b: BaseFeature | None) -> bool:
    """
    General function for determining whether feature 'a' and feature 'b' overlap.
    """
    if not (a and b):
        return False
    return intervals_overlap(a, b)


def feature_contains(a: BaseFeature | None, b: BaseFeature | None) -> bool:
    """
    General function for determining whether feature 'a' contains
    feature 'b'.  Note that both features must have 'minpos'
    and 'maxpos' attributes.
    """
    if not (a and b):
        return False
    return interval_contains(a, b)


def feature_search(
    features: Sequence[BaseFeature],
    query: BaseFeature,
    lo: int = 0,
    hi: int | None = None,
    overlap_window: int = 8,
) -> BaseFeature:
    """
    Bisect-based search through a sorted feature list.

    Returns either the feature that contains ``query`` or the feature that
    would immediately precede it.
    """
    if not features:
        raise ValueError("Cannot search an empty feature list")

    if hi is None:
        hi = len(features) - 1

    lo = max(0, lo)
    hi = min(hi, len(features) - 1)
    if lo > hi:
        raise ValueError("Invalid search bounds")

    index = InMemoryIntervalIndex(features)
    return index.predecessor_or_containing(
        query,
        lo=lo,
        hi=hi,
        overlap_window=overlap_window,
    )


def _resolve_gene_ptr(
    explicit_gene: Gene | None,
    parent: str | Gene | None,
    *,
    context: str,
) -> Gene:
    """Resolve parent gene for GTF serialization without requiring object back-links."""
    if explicit_gene is not None:
        return explicit_gene
    if isinstance(parent, Gene):
        return parent
    raise RuntimeError(context)


@dataclass(slots=True, eq=False)
class BaseFeature:
    feature_type: str | RecordType
    start_pos: InitVar[int]
    end_pos: InitVar[int]
    chromosome: str
    strand: str
    attr: InitVar[AttrMap | None] = None
    parent: str | Gene | None = field(init=False, default=None)
    minpos: int = field(init=False)
    maxpos: int = field(init=False)
    attributes: dict[str, str] = field(init=False, default_factory=dict)

    def __post_init__(
        self,
        start_pos: int,
        end_pos: int,
        attr: AttrMap | None,
    ) -> None:
        self.minpos = min(start_pos, end_pos)
        self.maxpos = max(start_pos, end_pos)
        self.attributes = (
            {} if attr is None else {str(key): str(value) for key, value in attr.items()}
        )

    def acceptor(self) -> int:
        """
        Returns the location where the acceptor dimer begins.

        Contract: use canonical splice-dimer assumptions with 1-based
        chromosome coordinates. On '+' this is exon-start minus 2 nt.
        On '-' this is the exact exon-start.
          + Example: AG|CGTATTC
          - Example: GAATACG|CT (reverse-complement)
        """
        return self.start() - SPLICE_DIMER_OFFSET if self.strand == "+" else self.start()

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
        )

    def contains(self, pos: int, strand: str) -> bool:
        return strand == self.strand and self.minpos <= pos <= self.maxpos

    def detail_string(self) -> str:
        return (
            f"Feature: {self.feature_type}\nChromosome: {self.chromosome}\n"
            f"Start: {self.start()}; End: {self.end()}; Strand: '{self.strand}'"
        )

    def donor(self) -> int:
        """
        Returns the location where the donor dimer begins.

        Contract: use canonical splice-dimer assumptions with 1-based
        chromosome coordinates. On '+' this is the exon end. On '-'
        this is exon-end minus 2 nt.
          + strand example: CGTATTC|GT
          - strand example: AC|GAATACG
        """
        return self.end() if self.strand == "+" else self.end() - SPLICE_DIMER_OFFSET

    def end(self) -> int:
        return self.maxpos if self.strand == "+" else self.minpos

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (
            self.minpos == other.minpos
            and self.maxpos == other.maxpos
            and self.strand == other.strand
            and self.chromosome == other.chromosome
        )

    def gff_string(self, alt_attributes: ExtraAttrMap | None = None) -> str:
        """Returns a GFF-formatted representation of the feature."""
        attrs = dict(self.attributes)
        if alt_attributes is not None:
            attrs.update({str(k): str(v) for k, v in alt_attributes.items()})

        return (
            f"{self.chromosome.capitalize()}\t{GFF_ID}\t{self.feature_type}\t"
            f"{self.minpos}\t{self.maxpos}\t.\t{self.strand}\t.\t{dict_to_gff(attrs)}"
        )

    def gtf_string(self, transcript: str, gene_ptr: Gene, exon_id: int) -> str:
        """Returns a GTF-formatted representation of the feature."""
        attrs = {str(k).lower(): str(v) for k, v in self.attributes.items()}
        attrs[GTF_GENE_ID] = gene_ptr.id
        attrs[GTF_GENE_NAME] = gene_ptr.name
        attrs[GTF_TRANSCRIPT] = transcript
        attrs[GTF_EXON_ID] = str(exon_id)
        try:
            source = attrs[GTF_SOURCE]
        except KeyError:
            source = GFF_ID
        return (
            f"{self.chromosome.capitalize()}\t{source}\t{self.feature_type}\t"
            f"{self.minpos}\t{self.maxpos}\t.\t{self.strand}\t.\t{dict_to_gtf(attrs)}"
        )

    def __hash__(self) -> int:
        return hash((self.minpos, self.maxpos, self.strand, self.chromosome))

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def set_parent(self, id: str | Gene) -> None:
        self.parent = id

    def start(self) -> int:
        return self.minpos if self.strand == "+" else self.maxpos

    def __str__(self) -> str:
        return f"{self.chromosome} {self.feature_type}: {self.start()}-{self.end()} ({self.strand})"


class TranscriptRegion(BaseFeature):
    def __init__(
        self,
        feature_type: str | RecordType,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        super().__init__(feature_type, start, end, chromosome, strand, attr)

    def add_parent(self, isoform: Isoform | Mrna) -> None:
        # Compatibility no-op: transcript regions no longer track child->parent object links.
        _ = isoform

    def __str__(self) -> str:
        return f"{self.feature_type} {self.start()}-{self.end()}({self.strand})"


class Exon(TranscriptRegion):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        super().__init__(RecordType.EXON, start, end, chromosome, strand, attr)


class Isoform(BaseFeature):
    def __init__(
        self,
        id: str,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, ISOFORM_TYPE, start, end, chromosome, strand, attr)
        self.id: str = id
        self.features: list[BaseFeature] = []
        self.exons: list[Exon] = []
        self.exon_map: dict[tuple[int, int], Exon] = {}

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, Isoform):
            return NotImplemented
        return self.feature_type == other.feature_type and self.id == other.id

    def acceptor_list(self) -> list[int]:
        """
        Returns a list of acceptor positions for the isoform.
        """
        ordered_exons = sorted(self.exons, key=feature_sort_key, reverse=(self.strand == "-"))
        return [e.acceptor() for e in ordered_exons[1:]]

    def add_exon(self, exon: Exon) -> bool:
        """
        Adds an exon to the isoform if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if exon.strand != self.strand:
            raise ValueError(
                f"Exon strand '{exon.strand}' does not match isoform "
                f"strand '{self.strand}' for {self.id}"
            )
        if exon.chromosome != self.chromosome:
            raise ValueError(
                f"Exon chromosome '{exon.chromosome}' does not match "
                f"isoform chromosome '{self.chromosome}' for {self.id}"
            )

        exon_tuple = (exon.minpos, exon.maxpos)
        if exon_tuple in self.exon_map:
            return False

        self.exons.append(exon)
        self.exon_map[exon_tuple] = exon
        self.minpos = min(self.minpos, exon.minpos)
        self.maxpos = max(self.maxpos, exon.maxpos)
        return True

    def add_feature(self, feature: BaseFeature) -> None:
        if feature.strand != self.strand:
            raise ValueError(
                f"ERROR: feature strand '{feature.strand}' does not match "
                f"form strand '{self.strand}'"
            )

        if feature.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: feature chromosome '{feature.chromosome}' does not "
                f"match form chromosome '{self.chromosome}'"
            )

        feature.set_parent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def detail_string(self) -> str:
        return (
            f"Isoform {self.id}\nStart: {self.start()}; End: {self.end()}; "
            f"Strand: {self.strand}\nExons: [{self.exon_string()}]\n"
        )

    def donor_list(self) -> list[int]:
        """
        Returns a list of donor positions for the isoform.
        """
        ordered_exons = sorted(self.exons, key=feature_sort_key, reverse=(self.strand == "-"))
        return [e.donor() for e in ordered_exons[:-1]]

    def exon_string(self) -> str:
        return ",".join([str(e) for e in self.exons])

    def get_feature_list(self, feature_type: str | RecordType) -> list[BaseFeature]:
        return [f for f in self.features if f.feature_type == feature_type]

    def gff_strings(self) -> list[str]:
        result: list[str] = []
        # Attributes depend on where data originated
        gff_attr = {PARENT_FIELD: self.id}
        gtf_attr = {
            PARENT_FIELD: self.id,
            GTF_TRANSCRIPT: self.id,
            GTF_TRANSNAME: self.id,
            GTF_PROTEIN_ID: self.id,
        }
        for exon in self.exons:
            if GTF_TRANSCRIPT in exon.attributes:
                result.append(exon.gff_string(alt_attributes=gtf_attr))
            else:
                result.append(exon.gff_string(alt_attributes=gff_attr))
        return result

    def gtf_strings(self, gene_ptr: Gene | None = None) -> list[str]:
        result: list[str] = []
        parent_gene = _resolve_gene_ptr(
            gene_ptr,
            self.parent,
            context="Isoform parent gene is required before writing GTF",
        )
        # Compatibility policy: GTF remains genomic-ascending on both strands.
        exon_list = sorted(self.exons, key=gtf_feature_sort_key)
        for i in range(len(exon_list)):
            exon = exon_list[i]
            result.append(exon.gtf_string(self.id, parent_gene, i + 1))
        return result

    def sorted_exons(self) -> list[Exon]:
        """
        Sorts the exons in an isoform based on its strand and
        returns the sorted list of exon objects.
        """
        return sorted(self.exons, key=feature_sort_key, reverse=(self.strand == "-"))

    def sorted_introns(self) -> list[tuple[int, int]]:
        """
        Returns a list of intron (donor,acceptor) tuples sorted based on strand.
        """
        result = []
        exons = self.sorted_exons()
        prev = exons[0]
        for e in exons[1:]:
            result.append((prev.donor(), e.acceptor()))
            prev = e
        return result

    def __str__(self) -> str:
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons)} exons/cds, "
            f"range {self.minpos} to {self.maxpos}"
        )


# Just as exons are part of an isoform, CDS/UTR regions are part of an Mrna sequence:
class CDS(TranscriptRegion):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        super().__init__(RecordType.CDS, start, end, chromosome, strand, attr)

    def __lt__(self, other: object) -> bool:
        """Compare CDS records, including feature type for tied intervals."""
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos, self.feature_type) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
            other.feature_type,
        )

    def __eq__(self, o: object) -> bool:
        """Compare CDS equality by type and genomic location."""
        if not isinstance(o, BaseFeature):
            return NotImplemented
        return (
            self.feature_type == o.feature_type
            and self.minpos == o.minpos
            and self.maxpos == o.maxpos
            and self.strand == o.strand
            and self.chromosome == o.chromosome
        )

    def __str__(self) -> str:
        return f"{self.feature_type} {self.start()}-{self.end()}({self.strand})"


# We treat UTR records the same way as CDS records
class FpUtr(TranscriptRegion):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        super().__init__(RecordType.FIVE_PRIME_UTR, start, end, chromosome, strand, attr)


class TpUtr(TranscriptRegion):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        super().__init__(RecordType.THREE_PRIME_UTR, start, end, chromosome, strand, attr)


class Mrna(Isoform):
    """
    An Mrna acts like an isoform in that it is associated with a parent gene
    and contains a number of coding sequences (CDS).
    """

    def __init__(
        self,
        id: str,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.MRNA, start, end, chromosome, strand, attr)
        self.id: str = id
        self.exons: list[Exon] = []
        self.features: list[BaseFeature] = []
        self.cds: list[TranscriptRegion] = []
        self.cds_map: dict[tuple[int, int], TranscriptRegion] = {}
        self.start_codon_pos: tuple[int, int] | None = None
        self.end_codon_pos: tuple[int, int] | None = None

    def acceptor_list(self) -> list[int]:
        """Returns a list of acceptor positions for the Mrna."""
        ordered_cds = sorted(self.cds, key=feature_sort_key, reverse=(self.strand == "-"))
        return [c.acceptor() for c in ordered_cds[1:]]

    def add_cds(self, cds: TranscriptRegion) -> bool:
        """
        Adds a CDS to the Mrna if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        if cds.strand != self.strand:
            raise ValueError(
                f"ERROR: CDS strand '{cds.strand}' does not match gene strand '{self.strand}'"
            )
        if cds.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: CDS chromosome '{cds.chromosome}' does not match "
                f"gene chromosome '{self.chromosome}'"
            )

        cds_tuple = (cds.minpos, cds.maxpos)
        if cds_tuple in self.cds_map:
            return False
        self.cds_map[cds_tuple] = cds

        self.cds.append(cds)
        return True

    def donor_list(self) -> list[int]:
        """Returns a list of donor positions for the Mrna."""
        ordered_cds = sorted(self.cds, key=feature_sort_key, reverse=(self.strand == "-"))
        return [c.donor() for c in ordered_cds[:-1]]

    def end_codon(self) -> tuple[int, int] | None:
        """Returns the end codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.end_codon_pos:
            self.find_codons()
        return self.end_codon_pos

    def find_codons(self) -> None:
        """Infers a transcript's start and end codon positions based on
        the relative positions of UTR and CDS records."""
        if self.end_codon_pos and self.start_codon_pos:
            return
        if not self.cds:
            return
        ordered_cds = sorted(self.cds, key=feature_sort_key, reverse=(self.strand == "-"))
        prev = ordered_cds[0]
        for c in ordered_cds[1:]:
            if (
                not self.start_codon_pos
                and prev.feature_type == RecordType.FIVE_PRIME_UTR
                and c.feature_type == RecordType.CDS
            ):
                self.start_codon_pos = (
                    (c.minpos, c.minpos + 2) if self.strand == "+" else (c.maxpos - 2, c.maxpos)
                )
            elif (
                not self.end_codon_pos
                and prev.feature_type == RecordType.CDS
                and c.feature_type == RecordType.THREE_PRIME_UTR
            ):
                self.end_codon_pos = (
                    (prev.maxpos - 2, prev.maxpos)
                    if self.strand == "+"
                    else (prev.minpos, prev.minpos + 2)
                )
            prev = c

    def infer_codons(self) -> None:
        """This method will infer start and end codons even when there are no
        UTR records for a transcript.  It is best to use this only after all
        data have been loaded for a gene."""
        # No CDS records --> no way to infer codons
        if not self.cds:
            return

        # First try UTR inference:
        self.find_codons()

        # If either codon is missing, assume the CDS endpoints
        # represent start/stop codons.
        ordered_cds = sorted(self.cds, key=feature_sort_key, reverse=(self.strand == "-"))
        if not self.start_codon_pos:
            c = ordered_cds[0]
            self.start_codon_pos = (
                (c.minpos, c.minpos + 2) if self.strand == "+" else (c.maxpos - 2, c.maxpos)
            )

        if not self.end_codon_pos:
            c = ordered_cds[-1]
            self.end_codon_pos = (
                (c.maxpos - 2, c.maxpos) if self.strand == "+" else (c.minpos, c.minpos + 2)
            )

    def get_utrs(self) -> list[TranscriptRegion]:
        """Returns a list of all UTR records in the Mrna object."""
        return [
            c
            for c in self.cds
            if c.feature_type in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]
        ]

    def gff_strings(self) -> list[str]:
        result = [self.gff_string()]
        gff_attr = {PARENT_FIELD: self.id}
        gtf_attr = {
            PARENT_FIELD: self.id,
            GTF_TRANSCRIPT: self.id,
            GTF_TRANSNAME: self.id,
            GTF_PROTEIN_ID: self.id,
        }
        for c in self.cds:
            if GTF_TRANSCRIPT in c.attributes:
                result.append(c.gff_string(alt_attributes=gtf_attr))
            else:
                result.append(c.gff_string(alt_attributes=gff_attr))
        return result

    def gtf_start_codon(self, gene_ptr: Gene | None = None) -> str | None:
        if self.start_codon_pos:
            parent_gene = _resolve_gene_ptr(
                gene_ptr,
                self.parent,
                context="Mrna parent gene is required before writing start_codon",
            )
            return (
                f"{self.chromosome}\t{GFF_ID}\tstart_codon\t{self.start_codon_pos[0]}"
                f"\t{self.start_codon_pos[1]}\t.\t{self.strand}\t.\tgene_id "
                f'"{parent_gene.id}"; transcript_id "{self.id}"'
            )
        return None

    def gtf_stop_codon(self, gene_ptr: Gene | None = None) -> str | None:
        if self.end_codon_pos:
            parent_gene = _resolve_gene_ptr(
                gene_ptr,
                self.parent,
                context="Mrna parent gene is required before writing stop_codon",
            )
            return (
                f"{self.chromosome}\t{GFF_ID}\tstop_codon\t{self.end_codon_pos[0]}"
                f"\t{self.end_codon_pos[1]}\t.\t{self.strand}\t.\tgene_id "
                f'"{parent_gene.id}"; transcript_id "{self.id}"'
            )
        return None

    def gtf_strings(self, gene_ptr: Gene | None = None) -> list[str]:
        """
        Returns GTF strings for all elements of the Mrna transcript, including
        start/stop codon locations.
        """
        result = []
        parent_gene = _resolve_gene_ptr(
            gene_ptr,
            self.parent,
            context="Mrna parent gene is required before writing GTF",
        )
        codon_string = (
            self.gtf_start_codon(parent_gene)
            if self.strand == "+"
            else self.gtf_stop_codon(parent_gene)
        )
        if codon_string:
            result.append(codon_string)

        # Compatibility policy: GTF remains genomic-ascending on both strands.
        cds_list = sorted(self.cds, key=gtf_feature_sort_key)
        for i in range(len(cds_list)):
            c = cds_list[i]
            result.append(c.gtf_string(self.id, parent_gene, i + 1))

        codon_string = (
            self.gtf_stop_codon(parent_gene)
            if self.strand == "+"
            else self.gtf_start_codon(parent_gene)
        )
        if codon_string:
            result.append(codon_string)

        return result

    def sorted_cds(self) -> list[TranscriptRegion]:
        """
        Sorts the CDS in an Mrna based on its strand and returns the sorted list of CDS objects.
        """
        return sorted(self.cds, key=feature_sort_key, reverse=(self.strand == "-"))

    def sorted_exons(self, *, minintron: int = 2) -> list[Exon]:
        """
        Infers an exon list from a CDS list and sorts the list based on strand.  Usually
        this means a 5' UTR record abuts a CDS record or a CDS record abuts a 3' UTR record,
        in which case an expanded exon represents both.
        """
        cds_list = self.sorted_cds()
        result: list[Exon] = []
        for i in range(len(cds_list)):
            cds = cds_list[i]
            if len(result) > 0 and abs(cds.start() - result[-1].end()) < minintron:
                prev = result[-1]
                result[-1] = Exon(prev.start(), cds.end(), self.chromosome, self.strand)
            else:
                result.append(Exon(cds.start(), cds.end(), self.chromosome, self.strand))
        # Revise feature types to indicate CDS instead of exon
        for e in result:
            e.feature_type = RecordType.CDS
        return result

    def start_codon(self) -> tuple[int, int] | None:
        """Returns the start codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.start_codon_pos:
            self.find_codons()
        return self.start_codon_pos

    def __str__(self) -> str:
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.cds)} exons/cds, "
            f"range {self.minpos} to {self.maxpos}"
        )


class Gene(BaseFeature):
    def __init__(
        self,
        id: str,
        note: str | None,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        name: str | None = None,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.GENE, start, end, chromosome, strand, attr)
        self.id: str = id
        self.name: str = name if name is not None else id
        self.note: str | None = note
        self.isoforms: dict[str, Isoform] = {}
        self.mrna: dict[str, Mrna] = {}
        self.exons: list[Exon] = []
        self.cds: list[TranscriptRegion] = []
        self.exon_map: dict[tuple[int, int], Exon] = {}
        self.cds_map: dict[tuple[str | RecordType, int, int], TranscriptRegion] = {}

        # Codons for all transcripts/Mrna, one entry per transcript
        self.start_codon_map: dict[str, tuple[int, int] | None] = {}
        self.end_codon_map: dict[str, tuple[int, int] | None] = {}

        # All other features associated with genes in an annotation file, such as:
        #    3'/5' UTRs, Mrna, miRNA, siRNA, tRNA, rRNA, ncRNA, snRNA, snoRNA
        self.features: list[BaseFeature] = []

    def acceptor_list(self) -> list[int]:
        """
        Returns a list of acceptors for this gene.
        """
        acceptor_set = set()
        for transcript in self._iter_transcripts():
            acceptor_set.update(transcript.acceptor_list())
        return sorted(acceptor_set, reverse=(self.strand == "-"))

    def _iter_isoforms(self) -> Iterable[Isoform]:
        return self.isoforms.values()

    def _iter_mrna_records(self) -> Iterable[Mrna]:
        return self.mrna.values()

    def _iter_transcripts(self) -> Iterable[Isoform | Mrna]:
        yield from self._iter_isoforms()
        yield from self._iter_mrna_records()

    def add_cds(self, new_mrna: Mrna, new_cds: TranscriptRegion) -> bool:
        """
        Adds a CDS to the gene if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        result = False
        cds_tuple = (new_cds.feature_type, new_cds.minpos, new_cds.maxpos)
        try:
            cds = self.cds_map[cds_tuple]
        except KeyError:
            cds = new_cds
            self.cds.append(cds)
            self.cds_map[cds_tuple] = cds
            self.minpos = min(self.minpos, cds.minpos)
            self.maxpos = max(self.maxpos, cds.maxpos)
            result = True

        mrna = self.add_mrna(new_mrna)
        mrna.add_cds(cds)
        self.start_codon_map[mrna.id] = mrna.start_codon()
        self.end_codon_map[mrna.id] = mrna.end_codon()
        return result

    def add_exon(self, new_isoform: Isoform, new_exon: Exon) -> bool:
        """
        Adds an exon to the gene if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if new_isoform is None:
            raise ValueError(f"Illegal null isoform in exon {new_exon}")
        result = False
        pos_tuple = (new_exon.minpos, new_exon.maxpos)
        try:
            exon = self.exon_map[pos_tuple]
        except KeyError:
            exon = new_exon
            self.exons.append(exon)
            self.exon_map[pos_tuple] = exon
            self.minpos = min(self.minpos, exon.minpos)
            self.maxpos = max(self.maxpos, exon.maxpos)
            result = True

        isoform = self.add_isoform(new_isoform)
        isoform.add_exon(exon)
        return result

    def add_feature(self, feature: BaseFeature) -> None:
        if feature.strand != self.strand:
            raise ValueError(
                f"ERROR: feature strand '{feature.strand}' does not match "
                f"gene strand '{self.strand}'"
            )

        if feature.chromosome != self.chromosome:
            raise ValueError(
                f"ERROR: feature chromosome '{feature.chromosome}' does not "
                f"match gene chromosome '{self.chromosome}'"
            )

        feature.set_parent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def add_isoform(self, isoform: Isoform) -> Isoform:
        isoform.set_parent(self.id)
        return self.isoforms.setdefault(isoform.id, isoform)

    def add_mrna(self, mrna: Mrna) -> Mrna:
        mrna.set_parent(self.id)
        return self.mrna.setdefault(mrna.id, mrna)

    def detail_string(self) -> str:
        result = (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons) + len(self.cds)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}\n"
        )
        result += "\n  ".join([x.detail_string() for x in self.isoforms.values()])
        return result

    def donor_list(self) -> list[int]:
        """
        Returns a list of donors for this gene.
        """
        donor_set = set()
        for transcript in self._iter_transcripts():
            donor_set.update(transcript.donor_list())
        return sorted(donor_set, reverse=(self.strand == "-"))

    def end_codons(self) -> dict[str, tuple[int, int] | None]:
        """Returns a dictionary of splice forms and their associated end codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.end_codon_map

    def get_feature_list(self, feature_type: str | RecordType) -> list[BaseFeature]:
        return [f for f in self.features if f.feature_type == feature_type]

    def get_introns(self) -> set[tuple[int, int]]:
        """
        Returns a list of duples containing start/end positions of introns in this gene.
        """
        result: set[tuple[int, int]] = set()
        for transcript in self._iter_transcripts():
            exons = transcript.sorted_exons()
            for i in range(1, len(exons)):
                result.add((exons[i - 1].end(), exons[i].start()))
        return result

    def get_isoform(self, id: str) -> Isoform:
        return self.isoforms[id]

    def get_junctions(self) -> set[tuple[int, int]]:
        """
        Returns a list of all known splice junctions for this gene
        based on its isoform list.  List contains only unique
        donor-acceptor duples.
        """
        result: set[tuple[int, int]] = set()
        for iid in self.isoforms:
            iso = self.isoforms[iid]
            exons = sorted(iso.exons, key=feature_sort_key, reverse=(self.strand == "-"))
            for i in range(1, len(exons)):
                result.add((exons[i - 1].donor(), exons[i].acceptor()))
        return result

    def gff_strings(self) -> str:
        """Returns a GFF string representation of the gene record plus
        all elements within the gene."""
        string_list = [self.gff_string(alt_attributes={ID_FIELD: self.id})]
        iso_set = set(self.isoforms.keys())
        mrna_set = set(self.mrna.keys())
        common_keys = iso_set & mrna_set
        iso_keys = iso_set - mrna_set
        mrna_keys = mrna_set - iso_set
        all_keys = sorted(iso_set | mrna_set)
        for k in all_keys:
            if k in common_keys:
                all_exons = self.isoforms[k].exons + self.mrna[k].cds
                all_exons.sort(key=feature_sort_key, reverse=(self.strand == "-"))
                string_list.append(self.mrna[k].gff_string())
                gff_attr = {PARENT_FIELD: k}
                gtf_attr = {
                    PARENT_FIELD: k,
                    GTF_TRANSCRIPT: k,
                    GTF_TRANSNAME: k,
                    GTF_PROTEIN_ID: k,
                }
                # string_list += [e.gff_string(alt_attributes={PARENT_FIELD:k}) for e in all_exons]
                for e in all_exons:
                    if GTF_TRANSCRIPT in e.attributes:
                        string_list.append(e.gff_string(alt_attributes=gtf_attr))
                    else:
                        string_list.append(e.gff_string(alt_attributes=gff_attr))

            elif k in iso_keys:
                string_list += self.isoforms[k].gff_strings()
            elif k in mrna_keys:
                string_list += self.mrna[k].gff_strings()

        return "\n".join(string_list)

    def gtf_strings(self) -> str:
        """Returns a GTF string representation of the gene record plus
        all elements within the gene."""
        string_list = []
        iso_set = set(self.isoforms.keys())
        mrna_set = set(self.mrna.keys())
        common_keys = iso_set & mrna_set
        iso_keys = iso_set - mrna_set
        mrna_keys = mrna_set - iso_set
        all_keys = sorted(iso_set | mrna_set)
        for k in all_keys:
            if k in common_keys:
                all_exons = self.isoforms[k].exons + self.mrna[k].cds
                all_exons.sort(key=gtf_feature_sort_key)

                # Ensure that there are codons to write
                self.mrna[k].infer_codons()

                codon_string = (
                    self.mrna[k].gtf_start_codon(self)
                    if self.strand == "+"
                    else self.mrna[k].gtf_stop_codon(self)
                )
                if codon_string:
                    string_list.append(codon_string)

                exon_counter = 0
                cds_counter = 0
                for i in range(len(all_exons)):
                    item = all_exons[i]
                    if item.feature_type == RecordType.EXON:
                        exon_counter += 1
                        string_list.append(item.gtf_string(k, self, exon_counter))
                    elif item.feature_type == RecordType.CDS:
                        cds_counter += 1
                        string_list.append(item.gtf_string(k, self, cds_counter))

                codon_string = (
                    self.mrna[k].gtf_stop_codon(self)
                    if self.strand == "+"
                    else self.mrna[k].gtf_start_codon(self)
                )
                if codon_string:
                    string_list.append(codon_string)

            elif k in iso_keys:
                string_list += self.isoforms[k].gtf_strings(self)
            elif k in mrna_keys:
                string_list += self.mrna[k].gtf_strings(self)
        return "\n".join(string_list)

    def is_single_exon(self) -> bool:
        return len(self.exons) == 1

    def sorted_exons(self) -> list[Exon]:
        """Returns a list of all exons inferred by the gene model, sorted 5' to 3'."""
        tmp_set = set()
        for isoform in self._iter_isoforms():
            tmp_set.update(isoform.sorted_exons())

        # Avoid returning duplicates:
        stored = {(e.minpos, e.maxpos) for e in tmp_set}
        for mrna_rec in self._iter_mrna_records():
            tmp_set.update(e for e in mrna_rec.sorted_exons() if (e.minpos, e.maxpos) not in stored)

        result = list(tmp_set)
        result.sort(key=feature_sort_key, reverse=(self.strand == "-"))
        return result

    def start_codons(self) -> dict[str, tuple[int, int] | None]:
        """Returns a dictionary of splice forms and their associated start codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.start_codon_map

    def __str__(self) -> str:
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.cds) + len(self.exons)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}"
        )


class PseudoGene(Gene):
    def __init__(
        self,
        id: str,
        note: str | None,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        name: str | None = None,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.PSEUDOGENE, start, end, chromosome, strand, attr)
        self.id = id
        self.name = name if name is not None else id
        self.note = note
        self.features = []
        self.exons = []
        self.cds = []
        self.mrna = {}
        self.isoforms = {}
        self.exon_map = {}
        self.cds_map = {}

        self.start_codon_map = {}
        self.end_codon_map = {}

    def detail_string(self) -> str:
        return self.__str__()
