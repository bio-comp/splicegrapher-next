"""
Stores gene annotation information from a GFF3 annotation
file and provides methods for searching on the data.
"""

from __future__ import annotations

# TRYING A NEW WAY TO IDENTIFY GENES:
import os
from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass
from typing import Protocol, TextIO, TypeVar
from urllib.parse import unquote

from SpliceGrapher.core.enums import AttrKey, RecordType, Strand
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)
from SpliceGrapher.formats.parsers.gene_model_gff import load_gene_model_records
from SpliceGrapher.formats.writers.gene_model import (
    write_gff as write_gene_model_gff,
)
from SpliceGrapher.formats.writers.gene_model import (
    write_gtf as write_gene_model_gtf,
)
from SpliceGrapher.shared.progress import ProgressIndicator

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
VALID_STRANDS = set(Strand)

FORM_DELIMITERS = [".", "-", "_", ","]

# Genes may overlap or even contain one another, so binary search
# is performed only at a coarse resolution to get within this number
# of genes.  Then a linear search is performed with a slight margin.
# This makes it approximately O(lg n) + C
GENE_SEARCH_RANGE = 16
GENE_SEARCH_MARGIN = 8

# GeneModel splice-site helpers assume canonical GT/AG-style 2nt dimer boundaries.
# This remains a compatibility contract for existing callers and fixtures.
SPLICE_DIMER_OFFSET = 2

# GTF serialization contract: emit records in genomic coordinate order for both
# strands (ascending minpos, then maxpos) to preserve historical output shape.
GTF_ORDER_POLICY = "genomic_ascending"

GeneFilter = Callable[["Gene"], bool]
GffRecordSource = str | list[str] | set[str] | tuple[str, ...]
AttrValue = str
AttrMap = Mapping[str, AttrValue]
ExtraAttrMap = Mapping[str, AttrValue] | Mapping[AttrKey, AttrValue]
AttrT = TypeVar("AttrT")
GeneLikeT = TypeVar("GeneLikeT", bound="GeneLike")


class GeneLike(Protocol):
    minpos: int
    maxpos: int
    strand: str

    def contains(self, pos: int, strand: str) -> bool: ...


@dataclass(slots=True)
class IntervalQuery:
    """Simple interval wrapper for index overlap queries."""

    minpos: int
    maxpos: int


def gene_type_filter(g: Gene) -> bool:
    """Convenience filter for getting only 'gene' records."""
    return g.featureType == RecordType.GENE


def defaultGeneFilter(g: Gene) -> bool:
    """Default function for filtering genes from a list."""
    return True


def dict_to_gff(d: Mapping[str, str]) -> str:
    """Returns a string representation of a dictionary based on the GFF3 annotation format."""
    return ";".join([f"{k}={v}" for k, v in sorted(d.items()) if k != "parent"])


def dict_to_gtf(d: Mapping[str, str]) -> str:
    """Returns a string representation of a dictionary based on the GTF annotation format."""
    return "; ".join([f'{k} "{v}"' for k, v in sorted(d.items()) if k != "parent"])


def cdsFactory(
    recType: RecordType,
    startPos: int,
    endPos: int,
    chrName: str,
    strand: str,
    attr: AttrMap | None = None,
) -> CDS:
    """Simple factory method for creating CDS-type records."""
    attr = {} if attr is None else {str(key): str(value) for key, value in attr.items()}
    if recType == RecordType.CDS:
        return CDS(startPos, endPos, chrName, strand, attr)
    elif recType == RecordType.FIVE_PRIME_UTR:
        return FP_UTR(startPos, endPos, chrName, strand, attr)
    elif recType == RecordType.THREE_PRIME_UTR:
        return TP_UTR(startPos, endPos, chrName, strand, attr)
    else:
        raise ValueError(f"Illegal CDS record type: {recType}")


class Chromosome:
    """Class that encapsulates a chromosome GFF record."""

    def __init__(self, start: int, end: int, name: str) -> None:
        self.minpos = start
        self.maxpos = end
        self.name = name

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def __str__(self) -> str:
        return f"{self.name}: {self.minpos}-{self.maxpos}"

    def contains(self, pos: int) -> bool:
        return self.minpos <= pos <= self.maxpos

    def end(self) -> int:
        return self.maxpos

    def gffString(self) -> str:
        # Example: Chr1    TAIR9   chromosome  1   30427671    .   .   .   ID=Chr1;Name=Chr1
        nameStr = self.name.capitalize()
        return (
            f"{nameStr}\t{GFF_ID}\tchromosome\t{self.start()}\t{self.end()}\t.\t.\t."
            f"\tID={nameStr};Name={nameStr}"
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


def featureCmp(a: BaseFeature, b: BaseFeature) -> int:
    """
    General comparison function for sorting any features that have 'minpos'
    and 'maxpos' attributes.
    """
    if a.minpos == b.minpos:
        return a.maxpos - b.maxpos
    else:
        return a.minpos - b.minpos


def featureSortKey(feature: BaseFeature) -> tuple[int, int]:
    """Sort features by genomic interval."""
    return (feature.minpos, feature.maxpos)


def geneSortKey(gene: Gene) -> tuple[int, int, str]:
    """Sort genes by interval, then id for deterministic ties."""
    return (gene.minpos, gene.maxpos, gene.id)


def gtf_feature_sort_key(feature: BaseFeature) -> tuple[int, int]:
    """Sort key for GTF emission using the genomic-ascending compatibility policy."""
    return featureSortKey(feature)


def featureOverlaps(a: BaseFeature | None, b: BaseFeature | None) -> bool:
    """
    General function for determining whether feature 'a' and feature 'b' overlap.
    """
    if not (a and b):
        return False
    return intervals_overlap(a, b)


def featureContains(a: BaseFeature | None, b: BaseFeature | None) -> bool:
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


class BaseFeature:
    def __init__(
        self,
        featureType: str | RecordType,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        self.chromosome: str = chromosome
        self.strand: str = strand
        self.parent: str | Gene | None = None
        self.minpos = min(start, end)
        self.maxpos = max(start, end)
        self.featureType: str | RecordType = featureType
        self.attributes: dict[str, str] = (
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

    def __lt__(self, other: BaseFeature) -> bool:
        return (self.chromosome, self.minpos, self.maxpos) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
        )

    def contains(self, pos: int, strand: str) -> bool:
        return strand == self.strand and self.minpos <= pos <= self.maxpos

    def detailString(self) -> str:
        return (
            f"Feature: {self.featureType}\nChromosome: {self.chromosome}\n"
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

    def __eq__(self, other) -> bool:
        return (
            self.minpos == other.minpos
            and self.maxpos == other.maxpos
            and self.strand == other.strand
            and self.chromosome == other.chromosome
        )

    def gffString(self, altAttributes: ExtraAttrMap | None = None) -> str:
        """Returns a GFF-formatted representation of the feature."""
        attrs = dict(self.attributes)
        if altAttributes is not None:
            attrs.update({str(k): str(v) for k, v in altAttributes.items()})

        return (
            f"{self.chromosome.capitalize()}\t{GFF_ID}\t{self.featureType}\t"
            f"{self.minpos}\t{self.maxpos}\t.\t{self.strand}\t.\t{dict_to_gff(attrs)}"
        )

    def gtfString(self, transcript: str, genePtr: Gene, exonId: int) -> str:
        """Returns a GTF-formatted representation of the feature."""
        attrs = {str(k).lower(): str(v) for k, v in self.attributes.items()}
        attrs[GTF_GENE_ID] = genePtr.id
        attrs[GTF_GENE_NAME] = genePtr.name
        attrs[GTF_TRANSCRIPT] = transcript
        attrs[GTF_EXON_ID] = str(exonId)
        try:
            source = attrs[GTF_SOURCE]
        except KeyError:
            source = GFF_ID
        return (
            f"{self.chromosome.capitalize()}\t{source}\t{self.featureType}\t"
            f"{self.minpos}\t{self.maxpos}\t.\t{self.strand}\t.\t{dict_to_gtf(attrs)}"
        )

    def __hash__(self) -> int:
        return hash((self.minpos, self.maxpos, self.strand, self.chromosome))

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def setParent(self, id: str | Gene) -> None:
        self.parent = id

    def start(self) -> int:
        return self.minpos if self.strand == "+" else self.maxpos

    def __str__(self) -> str:
        return f"{self.chromosome} {self.featureType}: {self.start()}-{self.end()} ({self.strand})"


class Exon(BaseFeature):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ):
        BaseFeature.__init__(self, RecordType.EXON, start, end, chromosome, strand, attr)
        self.parents: list[Isoform | mRNA] = []

    def __str__(self) -> str:
        if self.parents:
            return (
                f"{self.featureType} {self.start()}-{self.end()}({self.strand}) "
                f"(isoforms: {','.join([x.id for x in self.parents])})"
            )
        else:
            return f"{self.featureType} {self.start()}-{self.end()}({self.strand})"

    def addParent(self, isoform: Isoform | mRNA) -> None:
        if isoform not in self.parents:
            self.parents.append(isoform)


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
        self.exonMap: dict[tuple[int, int], Exon] = {}

    def __eq__(self, other) -> bool:
        return self.featureType == other.featureType and self.id == other.id

    def acceptorList(self) -> list[int]:
        """
        Returns a list of acceptor positions for the isoform.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [e.acceptor() for e in self.exons[1:]]

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

        exonTuple = (exon.minpos, exon.maxpos)
        if exonTuple in self.exonMap:
            return False

        self.exons.append(exon)
        self.exonMap[exonTuple] = exon
        if self not in exon.parents:
            exon.parents.append(self)
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

        feature.setParent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def addFeature(self, feature: BaseFeature) -> None:
        self.add_feature(feature)

    def detailString(self) -> str:
        return (
            f"Isoform {self.id}\nStart: {self.start()}; End: {self.end()}; "
            f"Strand: {self.strand}\nExons: [{self.exonString()}]\n"
        )

    def donorList(self) -> list[int]:
        """
        Returns a list of donor positions for the isoform.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [e.donor() for e in self.exons[:-1]]

    def exonString(self) -> str:
        return ",".join([str(e) for e in self.exons])

    def get_feature_list(self, featureType: str | RecordType) -> list[BaseFeature]:
        return [f for f in self.features if f.featureType == featureType]

    def gffStrings(self) -> list[str]:
        result: list[str] = []
        # Attributes depend on where data originated
        gffAttr = {PARENT_FIELD: self.id}
        gtfAttr = {
            PARENT_FIELD: self.id,
            GTF_TRANSCRIPT: self.id,
            GTF_TRANSNAME: self.id,
            GTF_PROTEIN_ID: self.id,
        }
        for exon in self.exons:
            if GTF_TRANSCRIPT in exon.attributes:
                result.append(exon.gffString(altAttributes=gtfAttr))
            else:
                result.append(exon.gffString(altAttributes=gffAttr))
        return result

    def gtfStrings(self) -> list[str]:
        result: list[str] = []
        if not isinstance(self.parent, Gene):
            raise RuntimeError("Isoform parent gene is required before writing GTF")
        # Compatibility policy: GTF remains genomic-ascending on both strands.
        exonList = sorted(self.exons, key=gtf_feature_sort_key)
        for i in range(len(exonList)):
            exon = exonList[i]
            result.append(exon.gtfString(self.id, self.parent, i + 1))
        return result

    def sortedExons(self) -> list[Exon]:
        """
        Sorts the exons in an isoform based on its strand and
        returns the sorted list of exon objects.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return self.exons

    def sortedIntrons(self) -> list[tuple[int, int]]:
        """
        Returns a list of intron (donor,acceptor) tuples sorted based on strand.
        """
        result = []
        exons = self.sortedExons()
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


# Just as exons are part of an isoform, so CDS elements are part of an mRNA sequence:
class CDS(Exon):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.CDS, start, end, chromosome, strand, attr)
        self.parents: list[Isoform | mRNA] = []

    def __lt__(self, other: BaseFeature) -> bool:
        """Compare CDS records, including feature type for tied intervals."""
        return (self.chromosome, self.minpos, self.maxpos, self.featureType) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
            other.featureType,
        )

    def __eq__(self, o) -> bool:
        """Compare CDS equality by type and genomic location."""
        return (
            self.featureType == o.featureType
            and self.minpos == o.minpos
            and self.maxpos == o.maxpos
            and self.strand == o.strand
            and self.chromosome == o.chromosome
        )

    def __str__(self) -> str:
        return (
            f"{self.featureType} {self.start()}-{self.end()} "
            f"(isoforms: {','.join([x.id for x in self.parents])})"
        )


# We treat UTR records the same way as CDS records
class FP_UTR(CDS):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.FIVE_PRIME_UTR, start, end, chromosome, strand, attr)
        self.parents: list[Isoform | mRNA] = []


class TP_UTR(CDS):
    def __init__(
        self,
        start: int,
        end: int,
        chromosome: str,
        strand: str,
        attr: AttrMap | None = None,
    ) -> None:
        BaseFeature.__init__(self, RecordType.THREE_PRIME_UTR, start, end, chromosome, strand, attr)
        self.parents: list[Isoform | mRNA] = []


class mRNA(Isoform):
    """
    An mRNA acts like an isoform in that it is associated with a parent gene
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
        self.cds: list[CDS] = []
        self.cdsMap: dict[tuple[int, int], CDS] = {}
        self.start_codon: tuple[int, int] | None = None
        self.end_codon: tuple[int, int] | None = None

    def acceptorList(self) -> list[int]:
        """Returns a list of acceptor positions for the mRNA."""
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [c.acceptor() for c in self.cds[1:]]

    def add_cds(self, cds: CDS) -> bool:
        """
        Adds a CDS to the mRNA if it's unique.
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

        cdsTuple = (cds.minpos, cds.maxpos)
        if cdsTuple in self.cdsMap:
            return False
        self.cdsMap[cdsTuple] = cds

        if self not in cds.parents:
            cds.parents.append(self)
        self.cds.append(cds)
        return True

    def donorList(self) -> list[int]:
        """Returns a list of donor positions for the mRNA."""
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [c.donor() for c in self.cds[:-1]]

    def endCodon(self) -> tuple[int, int] | None:
        """Returns the end codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.end_codon:
            self.findCodons()
        return self.end_codon

    def findCodons(self) -> None:
        """Infers a transcript's start and end codon positions based on
        the relative positions of UTR and CDS records."""
        if self.end_codon and self.start_codon:
            return
        if not self.cds:
            return
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        prev = self.cds[0]
        for c in self.cds[1:]:
            if (
                not self.start_codon
                and prev.featureType == RecordType.FIVE_PRIME_UTR
                and c.featureType == RecordType.CDS
            ):
                self.start_codon = (
                    (c.minpos, c.minpos + 2) if self.strand == "+" else (c.maxpos - 2, c.maxpos)
                )
            elif (
                not self.end_codon
                and prev.featureType == RecordType.CDS
                and c.featureType == RecordType.THREE_PRIME_UTR
            ):
                self.end_codon = (
                    (prev.maxpos - 2, prev.maxpos)
                    if self.strand == "+"
                    else (prev.minpos, prev.minpos + 2)
                )
            prev = c

    def inferCodons(self) -> None:
        """This method will infer start and end codons even when there are no
        UTR records for a transcript.  It is best to use this only after all
        data have been loaded for a gene."""
        # No CDS records --> no way to infer codons
        if not self.cds:
            return

        # First try UTR inference:
        self.findCodons()

        # If either codon is missing, assume the CDS endpoints
        # represent start/stop codons.
        if not self.start_codon:
            c = self.cds[0]
            self.start_codon = (
                (c.minpos, c.minpos + 2) if self.strand == "+" else (c.maxpos - 2, c.maxpos)
            )

        if not self.end_codon:
            c = self.cds[-1]
            self.end_codon = (
                (c.maxpos - 2, c.maxpos) if self.strand == "+" else (c.minpos, c.minpos + 2)
            )

    def getUTRs(self) -> list[CDS]:
        """Returns a list of all UTR records in the mRNA object."""
        return [
            c
            for c in self.cds
            if c.featureType in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]
        ]

    def gffStrings(self) -> list[str]:
        result = [self.gffString()]
        gffAttr = {PARENT_FIELD: self.id}
        gtfAttr = {
            PARENT_FIELD: self.id,
            GTF_TRANSCRIPT: self.id,
            GTF_TRANSNAME: self.id,
            GTF_PROTEIN_ID: self.id,
        }
        for c in self.cds:
            if GTF_TRANSCRIPT in c.attributes:
                result.append(c.gffString(altAttributes=gtfAttr))
            else:
                result.append(c.gffString(altAttributes=gffAttr))
        return result

    def gtfStartCodon(self) -> str | None:
        if self.start_codon:
            if not isinstance(self.parent, Gene):
                raise RuntimeError("mRNA parent gene is required before writing start_codon")
            return (
                f"{self.chromosome}\t{GFF_ID}\tstart_codon\t{self.start_codon[0]}"
                f"\t{self.start_codon[1]}\t.\t{self.strand}\t.\tgene_id "
                f'"{self.parent.id}"; transcript_id "{self.id}"'
            )
        return None

    def gtfStopCodon(self) -> str | None:
        if self.end_codon:
            if not isinstance(self.parent, Gene):
                raise RuntimeError("mRNA parent gene is required before writing stop_codon")
            return (
                f"{self.chromosome}\t{GFF_ID}\tstop_codon\t{self.end_codon[0]}"
                f"\t{self.end_codon[1]}\t.\t{self.strand}\t.\tgene_id "
                f'"{self.parent.id}"; transcript_id "{self.id}"'
            )
        return None

    def gtfStrings(self) -> list[str]:
        """
        Returns GTF strings for all elements of the mRNA transcript, including
        start/stop codon locations.
        """
        result = []
        codonString = self.gtfStartCodon() if self.strand == "+" else self.gtfStopCodon()
        if codonString:
            result.append(codonString)
        if not isinstance(self.parent, Gene):
            raise RuntimeError("mRNA parent gene is required before writing GTF")

        # Compatibility policy: GTF remains genomic-ascending on both strands.
        cdsList = sorted(self.cds, key=gtf_feature_sort_key)
        for i in range(len(cdsList)):
            c = cdsList[i]
            result.append(c.gtfString(self.id, self.parent, i + 1))

        codonString = self.gtfStopCodon() if self.strand == "+" else self.gtfStartCodon()
        if codonString:
            result.append(codonString)

        return result

    def sortedCDS(self) -> list[CDS]:
        """
        Sorts the CDS in an mRNA based on its strand and returns the sorted list of CDS objects.
        """
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return self.cds

    def sortedExons(self, *, minintron: int = 2) -> list[Exon]:
        """
        Infers an exon list from a CDS list and sorts the list based on strand.  Usually
        this means a 5' UTR record abuts a CDS record or a CDS record abuts a 3' UTR record,
        in which case an expanded exon represents both.
        """
        cdsList = self.sortedCDS()
        result: list[Exon] = []
        for i in range(len(cdsList)):
            cds = cdsList[i]
            if len(result) > 0 and abs(cds.start() - result[-1].end()) < minintron:
                prev = result[-1]
                result[-1] = Exon(prev.start(), cds.end(), self.chromosome, self.strand)
            else:
                result.append(Exon(cds.start(), cds.end(), self.chromosome, self.strand))
        # Revise feature types to indicate CDS instead of exon
        for e in result:
            e.featureType = RecordType.CDS
        return result

    def startCodon(self) -> tuple[int, int] | None:
        """Returns the start codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.start_codon:
            self.findCodons()
        return self.start_codon

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
        self.mrna: dict[str, mRNA] = {}
        self.exons: list[Exon] = []
        self.cds: list[CDS] = []
        self.exonMap: dict[tuple[int, int], Exon] = {}
        self.cdsMap: dict[tuple[str | RecordType, int, int], CDS] = {}

        # Codons for all transcripts/mRNA, one entry per transcript
        self.start_codons: dict[str, tuple[int, int] | None] = {}
        self.end_codons: dict[str, tuple[int, int] | None] = {}

        # All other features associated with genes in an annotation file, such as:
        #    3'/5' UTRs, mRNA, miRNA, siRNA, tRNA, rRNA, ncRNA, snRNA, snoRNA
        self.features: list[BaseFeature] = []

    def acceptorList(self) -> list[int]:
        """
        Returns a list of acceptors for this gene.
        """
        acceptorSet = set()
        for transcript in self._iter_transcripts():
            acceptorSet.update(transcript.acceptorList())
        return sorted(list(acceptorSet), reverse=(self.strand == "-"))

    def _iter_isoforms(self) -> Iterable[Isoform]:
        return self.isoforms.values()

    def _iter_mrna_records(self) -> Iterable[mRNA]:
        return self.mrna.values()

    def _iter_transcripts(self) -> Iterable[Isoform | mRNA]:
        yield from self._iter_isoforms()
        yield from self._iter_mrna_records()

    def add_cds(self, new_mrna: mRNA, new_cds: CDS) -> bool:
        """
        Adds a CDS to the gene if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        result = False
        cdsTuple = (new_cds.featureType, new_cds.minpos, new_cds.maxpos)
        try:
            cds = self.cdsMap[cdsTuple]
        except KeyError:
            cds = new_cds
            self.cds.append(cds)
            self.cdsMap[cdsTuple] = cds
            self.minpos = min(self.minpos, cds.minpos)
            self.maxpos = max(self.maxpos, cds.maxpos)
            result = True

        mrna = self.add_mrna(new_mrna)
        mrna.add_cds(cds)
        self.start_codons[mrna.id] = mrna.startCodon()
        self.end_codons[mrna.id] = mrna.endCodon()
        return result

    def add_exon(self, new_isoform: Isoform, new_exon: Exon) -> bool:
        """
        Adds an exon to the gene if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if new_isoform is None:
            raise ValueError(f"Illegal null isoform in exon {new_exon}")
        result = False
        posTuple = (new_exon.minpos, new_exon.maxpos)
        try:
            exon = self.exonMap[posTuple]
        except KeyError:
            exon = new_exon
            self.exons.append(exon)
            self.exonMap[posTuple] = exon
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

        feature.setParent(self.id)
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)
        self.features.append(feature)

    def add_isoform(self, isoform: Isoform) -> Isoform:
        isoform.setParent(self)
        return self.isoforms.setdefault(isoform.id, isoform)

    def add_mrna(self, mrna: mRNA) -> mRNA:
        mrna.parent = self
        return self.mrna.setdefault(mrna.id, mrna)

    def detailString(self) -> str:
        result = (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons) + len(self.cds)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}\n"
        )
        result += "\n  ".join([x.detailString() for x in self.isoforms.values()])
        return result

    def donorList(self) -> list[int]:
        """
        Returns a list of donors for this gene.
        """
        donorSet = set()
        for transcript in self._iter_transcripts():
            donorSet.update(transcript.donorList())
        return sorted(list(donorSet), reverse=(self.strand == "-"))

    def endCodons(self) -> dict[str, tuple[int, int] | None]:
        """Returns a dictionary of splice forms and their associated end codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.end_codons

    def get_feature_list(self, featureType: str | RecordType) -> list[BaseFeature]:
        return [f for f in self.features if f.featureType == featureType]

    def getIntrons(self):
        """
        Returns a list of duples containing start/end positions of introns in this gene.
        """
        result = {}
        for transcript in self._iter_transcripts():
            exons = transcript.sortedExons()
            for i in range(1, len(exons)):
                key = (exons[i - 1].end(), exons[i].start())
                result[key] = 1

        return result.keys()

    def getIsoform(self, id: str) -> Isoform:
        return self.isoforms[id]

    def getJunctions(self):
        """
        Returns a list of all known splice junctions for this gene
        based on its isoform list.  List contains only unique
        donor-acceptor duples.
        """
        result = {}
        for iid in self.isoforms.keys():
            iso = self.isoforms[iid]
            exons = sorted(iso.exons, key=featureSortKey, reverse=(self.strand == "-"))
            for i in range(1, len(exons)):
                duple = (exons[i - 1].donor(), exons[i].acceptor())
                result[duple] = 1
        return result.keys()

    def gffStrings(self) -> str:
        """Returns a GFF string representation of the gene record plus
        all elements within the gene."""
        stringList = [self.gffString(altAttributes={ID_FIELD: self.id})]
        isoSet = set(self.isoforms.keys())
        mrnaSet = set(self.mrna.keys())
        commonKeys = isoSet & mrnaSet
        isoKeys = isoSet - mrnaSet
        mrnaKeys = mrnaSet - isoSet
        allKeys = sorted(isoSet | mrnaSet)
        for k in allKeys:
            if k in commonKeys:
                allExons = self.isoforms[k].exons + self.mrna[k].cds
                allExons.sort(key=featureSortKey, reverse=(self.strand == "-"))
                stringList.append(self.mrna[k].gffString())
                gffAttr = {PARENT_FIELD: k}
                gtfAttr = {
                    PARENT_FIELD: k,
                    GTF_TRANSCRIPT: k,
                    GTF_TRANSNAME: k,
                    GTF_PROTEIN_ID: k,
                }
                # stringList += [e.gffString(altAttributes={PARENT_FIELD:k}) for e in allExons]
                for e in allExons:
                    if GTF_TRANSCRIPT in e.attributes:
                        stringList.append(e.gffString(altAttributes=gtfAttr))
                    else:
                        stringList.append(e.gffString(altAttributes=gffAttr))

            elif k in isoKeys:
                stringList += self.isoforms[k].gffStrings()
            elif k in mrnaKeys:
                stringList += self.mrna[k].gffStrings()

        return "\n".join(stringList)

    def gtfStrings(self) -> str:
        """Returns a GTF string representation of the gene record plus
        all elements within the gene."""
        stringList = []
        isoSet = set(self.isoforms.keys())
        mrnaSet = set(self.mrna.keys())
        commonKeys = isoSet & mrnaSet
        isoKeys = isoSet - mrnaSet
        mrnaKeys = mrnaSet - isoSet
        allKeys = sorted(isoSet | mrnaSet)
        for k in allKeys:
            if k in commonKeys:
                allExons = self.isoforms[k].exons + self.mrna[k].cds
                allExons.sort(key=gtf_feature_sort_key)

                # Ensure that there are codons to write
                self.mrna[k].inferCodons()

                codonString = (
                    self.mrna[k].gtfStartCodon()
                    if self.strand == "+"
                    else self.mrna[k].gtfStopCodon()
                )
                if codonString:
                    stringList.append(codonString)

                eCtr = 0
                cCtr = 0
                for i in range(len(allExons)):
                    item = allExons[i]
                    if item.featureType == RecordType.EXON:
                        eCtr += 1
                        stringList.append(item.gtfString(k, self, eCtr))
                    elif item.featureType == RecordType.CDS:
                        cCtr += 1
                        stringList.append(item.gtfString(k, self, cCtr))

                codonString = (
                    self.mrna[k].gtfStopCodon()
                    if self.strand == "+"
                    else self.mrna[k].gtfStartCodon()
                )
                if codonString:
                    stringList.append(codonString)

            elif k in isoKeys:
                stringList += self.isoforms[k].gtfStrings()
            elif k in mrnaKeys:
                stringList += self.mrna[k].gtfStrings()
        return "\n".join(stringList)

    def isSingleExon(self) -> bool:
        return len(self.exons) == 1

    def sortedExons(self) -> list[Exon]:
        """Returns a list of all exons inferred by the gene model, sorted 5' to 3'."""
        tmpset = set()
        for isoform in self._iter_isoforms():
            tmpset.update(isoform.sortedExons())

        # Avoid returning duplicates:
        stored = set([(e.minpos, e.maxpos) for e in tmpset])
        for mrna_rec in self._iter_mrna_records():
            tmpset.update([e for e in mrna_rec.sortedExons() if (e.minpos, e.maxpos) not in stored])

        result = list(tmpset)
        result.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return result

    def startCodons(self) -> dict[str, tuple[int, int] | None]:
        """Returns a dictionary of splice forms and their associated start codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.start_codons

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
        self.exonMap = {}
        self.cdsMap = {}

        self.start_codons = {}
        self.end_codons = {}

    def detailString(self) -> str:
        return self.__str__()


class GeneModel:
    def __init__(
        self,
        gff_path: GffRecordSource | None,
        *,
        require_notes: bool = False,
        chromosomes: Sequence[str] | str | None = None,
        verbose: bool = False,
        ignore_errors: bool = False,
    ) -> None:
        """Instantiates a GeneModel object.  If a path is provided, this will
        load gene models from the given file."""
        # 2-dimensional model indexed by chromosome then gene
        self.allGenes: dict[str, Gene] = {}
        self.allChr: dict[str, Chromosome] = {}
        self.foundTypes: dict[RecordType, bool] = {}
        self.model: dict[str, dict[str, Gene]] = {}
        self.mRNAforms: dict[str, dict[str, mRNA]] = {}
        self.sorted: dict[str, dict[str, list[Gene]]] = {}
        self.interval_index: dict[str, dict[str, InMemoryIntervalIndex[Gene]]] = {}

        if gff_path:  # Load gene models from a file
            if isinstance(gff_path, str) and not os.path.exists(gff_path):
                raise ValueError(f"Gene model file not found: {gff_path}")

            self.load_gene_model(
                gff_path,
                require_notes=require_notes,
                chromosomes=chromosomes,
                verbose=verbose,
                ignore_errors=ignore_errors,
            )

            if not self.model:
                raise ValueError(f"No gene models found in {gff_path}")
            self.make_sorted_model()

    def __contains__(self, gene: str | Gene) -> bool:
        """Returns true if a gene is in the model; false otherwise."""
        return str(gene) in self.allGenes

    def add_chromosome(self, start: int, end: int, name: str) -> None:
        """Adds a chromosome to a gene model or updates the end points
        if the record already exists."""
        key = name.lower()
        try:
            rec = self.allChr[key]
            rec.minpos = min(rec.minpos, start)
            rec.maxpos = max(rec.maxpos, end)
        except KeyError:
            self.allChr[key] = Chromosome(start, end, name)
            self.model.setdefault(key, {})

    def add_gene(self, gene: Gene) -> None:
        """Adds a gene to a gene model.  Raises a ValueError if the
        gene has already been added."""
        if gene.id in self.model[gene.chromosome] or str(gene) in self.allGenes:
            raise ValueError(f"Gene {gene.id} already stored in gene model")

        self.model[gene.chromosome][gene.id] = gene
        self.allGenes[str(gene)] = gene

    def binary_search_genes(
        self,
        geneList: list[GeneLikeT],
        loc: int,
        lo: int,
        hi: int,
        pfx: str = "",
    ) -> tuple[GeneLikeT | None, GeneLikeT | None]:
        """Quasi-binary search through a gene list.  Performs binary search to a point,
        then performs linear search within a refined gene range.  Real binary search may
        not be possible since gene models may overlap or contain one within another."""
        _ = pfx
        if not geneList:
            return (None, None)

        lo_idx = lo
        hi_idx = hi
        while lo_idx <= hi_idx:
            if (hi_idx - lo_idx) <= GENE_SEARCH_RANGE:
                lo_bound = max(0, lo_idx - GENE_SEARCH_MARGIN)
                hi_bound = min(hi_idx + GENE_SEARCH_MARGIN, len(geneList) - 1)
                for idx in range(lo_bound, hi_bound + 1):
                    gene = geneList[idx]
                    if gene.contains(loc, gene.strand):
                        return (gene, gene)
                return (None, None)

            midpt = (hi_idx + lo_idx) // 2
            midGene = geneList[midpt]
            if midGene.contains(loc, midGene.strand):
                return (midGene, midGene)
            if loc > midGene.maxpos:
                lo_idx = midpt + 1
            elif loc < midGene.minpos:
                hi_idx = midpt - 1
            else:
                break
        return (None, None)

    def clean_name(self, s: str) -> str:
        """
        Some feature names include URL characters that we may wish to fix.
        """
        revised = unquote(s)
        revised = revised.replace(",", "")
        return revised.replace(" ", "-")

    def get_all_acceptors(
        self, geneFilter: GeneFilter = defaultGeneFilter
    ) -> dict[str, dict[str, set[int]]]:
        """
        Returns a dictionary of all known acceptor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model.keys():
            result[chrom] = self.get_known_acceptors(chrom, geneFilter)
        return result

    def get_all_donors(
        self, geneFilter: GeneFilter = defaultGeneFilter
    ) -> dict[str, dict[str, set[int]]]:
        """
        Returns a dictionary of all known donor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model.keys():
            result[chrom] = self.get_known_donors(chrom, geneFilter)
        return result

    def get_all_gene_ids(self, geneFilter: GeneFilter = defaultGeneFilter) -> list[str]:
        """Returns a list of ids for all genes stored."""
        return [g.id for g in self.allGenes.values() if geneFilter(g)]

    def get_all_genes(
        self,
        geneFilter: GeneFilter = defaultGeneFilter,
        *,
        verbose: bool = False,
    ) -> list[Gene]:
        """Returns a list of all genes stored."""
        indicator = ProgressIndicator(10000, verbose=verbose)
        result = []
        for g in self.allGenes.values():
            indicator.update()
            if geneFilter(g):
                result.append(g)
        indicator.finish()
        return result

    def get_annotation(
        self, key: str, annotDict: dict[str, str], default: AttrT | None = None
    ) -> str | AttrT | None:
        """
        Convenience method for retrieving a value from an annotation dictionary
        """
        try:
            return annotDict[key]
        except KeyError:
            return default

    def get_annotation_dict(self, s: str) -> dict[str, str]:
        """
        Parses a ';'-separated annotation string containing key-value pairs
        and returns them as a dictionary.
        """
        valStr = s.replace(" ", "")
        parts = valStr.split(";")
        result = {}
        for p in parts:
            if "=" not in p:
                continue
            key, value = p.split("=", 1)
            if not key:
                continue
            result[key] = value
        return result

    def get_chromosome(self, chrName: str) -> Chromosome | None:
        """Returns a simple record with basic chromosome information."""
        try:
            return self.allChr[chrName]
        except KeyError:
            return None

    def get_chromosomes(self) -> Iterable[str]:
        """Returns a list of all chromosomes represented in the model."""
        return self.model.keys()

    def get_feature_list(self, featureType: str | RecordType) -> list[BaseFeature]:
        """
        Returns a list of all features of the given type found in all genes.
        """
        result = []
        for gene in self.allGenes.values():
            fList = gene.get_feature_list(featureType)
            if fList:
                result += fList
        return result

    def get_gene(self, chrom: str, geneId: str) -> Gene | None:
        """
        Returns a gene from within a chromosome.  The gene will contain
        information on all exons within it.
        """
        try:
            return self.model[chrom.lower()][geneId]
        except KeyError:
            return None

    def get_gene_by_name(self, id: str) -> Gene | None:
        """
        Returns a gene with the given id if it exists.
        """
        for k in self.model.keys():
            try:
                return self.model[k][id.upper()]
            except KeyError:
                pass
        return None

    def get_gene_from_locations(
        self, chrom: str, startPos: int, endPos: int, strand: str
    ) -> Gene | None:
        """
        Finds the gene within the given chromosome that contains the given start
        and end positions.

        Uses the in-memory interval index built by ``make_sorted_model``.
        """
        chrom_key = chrom.lower()
        try:
            strand_index = self.interval_index[chrom_key][strand]
        except KeyError:
            raise KeyError(f"Key {chrom_key} not found in {','.join(self.sorted.keys())}")

        query = IntervalQuery(min(startPos, endPos), max(startPos, endPos))
        for gene in strand_index.overlaps(query):
            if gene.contains(startPos, strand) or gene.contains(endPos, strand):
                return gene
        return None

    def get_gene_records(
        self,
        chrom: str,
        geneFilter: GeneFilter = defaultGeneFilter,
        *,
        verbose: bool = False,
    ) -> list[Gene]:
        """
        Returns a list of all gene instances represented within a given chromosome.
        The gene list may be filtered by changing the geneFilter function.
        """
        try:
            # return [g for g in self.model[chrom.lower()].values() if geneFilter(g)]
            indicator = ProgressIndicator(10000, verbose=verbose)
            result = []
            for g in self.model[chrom.lower()].values():
                indicator.update()
                if geneFilter(g):
                    result.append(g)
            indicator.finish()
            return result
        except KeyError:
            return []

    def get_genes(self, chrom: str) -> list[str]:
        """
        Returns a list of all genes represented within a given chromosome.
        """
        try:
            return list(self.model[chrom.lower()].keys())
        except KeyError:
            return []

    def get_genes_in_range(
        self, chrom: str, minpos: int, maxpos: int, strand: str | None = None
    ) -> list[Gene]:
        """
        Returns a list of all gene instances represented within a given chromosome
        that overlap the given range.  If no strand is specified, this will
        return all genes on both strands that overlap the range.
        """
        result = []
        for g in self.get_gene_records(chrom):
            if g.maxpos < minpos or g.minpos > maxpos:
                continue
            if strand and g.strand != strand:
                continue
            result.append(g)
        return result

    def get_known_acceptors(
        self, chrom: str, geneFilter: GeneFilter = defaultGeneFilter
    ) -> dict[str, set[int]]:
        """
        Returns a dictionary of all known acceptor sites for the chromosome,
        indexed by strand.
        """
        result: dict[str, set[int]] = {"-": set(), "+": set()}
        for g in self.get_gene_records(chrom, geneFilter):
            result[g.strand].update(g.acceptorList())
        return result

    def get_known_donors(
        self, chrom: str, geneFilter: GeneFilter = defaultGeneFilter
    ) -> dict[str, set[int]]:
        """
        Returns a dictionary of all known donor sites for the chromosome,
        indexed by strand.
        """
        result: dict[str, set[int]] = {"-": set(), "+": set()}
        for g in self.get_gene_records(chrom, geneFilter):
            result[g.strand].update(g.donorList())
        return result

    def get_parent(
        self,
        s: str,
        chrom: str,
        searchGenes: bool = True,
        searchmRNA: bool = True,
    ) -> Gene | mRNA | None:
        """
        Parent identifiers are not stored in a consistent manner.  We may have
        'AT1G01160', 'AT1G01160.1' or '12345.AT1G01160' or possibly something else.
        This method looks for the most specific candidate name to identify a record's parent.
        """
        parent_key = s.upper()
        chrom_key = chrom.lower()

        if searchmRNA:
            chrom_forms = self.mRNAforms.get(chrom_key)
            if chrom_forms is not None and parent_key in chrom_forms:
                return chrom_forms[parent_key]

        if searchGenes:
            chrom_genes = self.model.get(chrom_key)
            if chrom_genes is not None and parent_key in chrom_genes:
                return chrom_genes[parent_key]
        return None

    def get_record_types(self) -> list[RecordType]:
        """Returns a list of all record types found in the input file."""
        return [k for k in self.foundTypes if self.foundTypes[k]]

    def isoform_dict(
        self,
        geneFilter: GeneFilter = defaultGeneFilter,
        *,
        verbose: bool = False,
    ) -> dict[str, set[str]]:
        """Returns a dictionary that maps gene names to their corresponding
        isoform identifiers.  Each gene is associated with a set of isoform ids."""
        result: dict[str, set[str]] = {}
        for g in self.get_all_genes(geneFilter=geneFilter, verbose=verbose):
            result[g.id] = set(list(g.mrna.keys()) + list(g.isoforms.keys()))
        return result

    def load_gene_model(
        self,
        gff_records: GffRecordSource,
        *,
        require_notes: bool = False,
        chromosomes: Sequence[str] | str | None = None,
        verbose: bool = False,
        ignore_errors: bool = False,
    ) -> None:
        """
        Reads a tab-delimited gene annotation GFF file and stores information
        on chromosomes, the genes within each chromosome and exons within each gene.

        Parameters:
          'gff_records'  - source of GFF records; may be a file path, a file stream
                           or a list/set of strings
          'require_notes' - require annotations for all gene records (default=False)
          'chromosomes'  - chromosome name or list of chromosomes to store (default=all)
          'verbose'      - provide verbose feedback (default=False)
          'ignore_errors' - ignore error conditions (default=False)
        """
        load_gene_model_records(
            self,
            gff_records,
            require_notes=require_notes,
            chromosomes=chromosomes,
            verbose=verbose,
            ignore_errors=ignore_errors,
        )

    def make_sorted_model(self) -> None:
        self.sorted = {}
        self.interval_index = {}
        for chrom in self.model.keys():
            self.sorted[chrom] = {}
            self.interval_index[chrom] = {}
            # Get a list of gene objects sorted by position
            geneList = sorted(
                self.model[chrom].values(),
                key=geneSortKey,
            )
            # Note: we accept unknown ('.') strands, but it's
            # up to calling routines to ask for them specifically
            self.sorted[chrom] = {"-": [], "+": [], ".": []}
            for g in geneList:
                self.sorted[chrom][g.strand].append(g)
            for strand, strand_genes in self.sorted[chrom].items():
                if strand_genes:
                    self.interval_index[chrom][strand] = InMemoryIntervalIndex(strand_genes)

    def write_gff(
        self,
        gff_path: str | TextIO,
        *,
        geneFilter: GeneFilter = defaultGeneFilter,
        geneSet: set[str] | list[str] | tuple[str, ...] | None = None,
        verbose: bool = False,
    ) -> None:
        """Writes a complete gene model out to a GFF file."""
        write_gene_model_gff(
            self,
            gff_path,
            gene_filter=geneFilter,
            gene_set=geneSet,
            verbose=verbose,
        )

    def write_gtf(
        self,
        gtf_path: str | TextIO,
        *,
        geneFilter: GeneFilter = defaultGeneFilter,
        verbose: bool = False,
    ) -> None:
        """Writes a complete gene model out to a GTF file."""
        write_gene_model_gtf(
            self,
            gtf_path,
            gene_filter=geneFilter,
            verbose=verbose,
        )
