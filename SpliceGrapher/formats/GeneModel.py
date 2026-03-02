"""
Stores gene annotation information from a GFF3 annotation
file and provides methods for searching on the data.
"""

from __future__ import annotations

# TRYING A NEW WAY TO IDENTIFY GENES:
import os
import sys
from urllib.parse import unquote

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AttrKey, RecordType, Strand
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.format_utils import comma_format
from SpliceGrapher.shared.process_utils import getAttribute
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


def gene_type_filter(g):
    """Convenience filter for getting only 'gene' records."""
    return g.featureType == RecordType.GENE


def defaultGeneFilter(g):
    """Default function for filtering genes from a list."""
    return True


def dictToGFF(d):
    """Returns a string representation of a dictionary based on the GFF3 annotation format."""
    return ";".join(["{}={}".format(k, v) for k, v in sorted(d.items()) if k != "parent"])


def dictToGTF(d):
    """Returns a string representation of a dictionary based on the GTF annotation format."""
    return "; ".join(['{} "{}"'.format(k, v) for k, v in sorted(d.items()) if k != "parent"])


def cdsFactory(recType, startPos, endPos, chrName, strand, attr=None):
    """Simple factory method for creating CDS-type records."""
    attr = {} if attr is None else attr
    if recType == RecordType.CDS:
        return CDS(startPos, endPos, chrName, strand, attr)
    elif recType == RecordType.FIVE_PRIME_UTR:
        return FP_UTR(startPos, endPos, chrName, strand, attr)
    elif recType == RecordType.THREE_PRIME_UTR:
        return TP_UTR(startPos, endPos, chrName, strand, attr)
    else:
        raise ValueError("Illegal CDS record type: {}".format(recType))


class Chromosome(object):
    """Class that encapsulates a chromosome GFF record."""

    def __init__(self, start, end, name):
        self.minpos = start
        self.maxpos = end
        self.name = name

    def __len__(self):
        return self.maxpos - self.minpos + 1

    def __str__(self):
        return f"{self.name}: {self.minpos}-{self.maxpos}"

    def contains(self, pos):
        return self.minpos <= pos <= self.maxpos

    def end(self):
        return self.maxpos

    def gffString(self):
        # Example: Chr1    TAIR9   chromosome  1   30427671    .   .   .   ID=Chr1;Name=Chr1
        nameStr = self.name.capitalize()
        return (
            f"{nameStr}\t{GFF_ID}\tchromosome\t{self.start()}\t{self.end()}\t.\t.\t."
            f"\tID={nameStr};Name={nameStr}"
        )

    def start(self):
        return self.minpos

    def update(self, feature):
        """
        Many species do not have 'chromosome' entries in their annotations,
        so we must infer chromosome boundaries from observed features.
        """
        if feature.chromosome.lower() != self.name.lower():
            raise ValueError(
                "Cannot use feature from {} to update {}".format(feature.chromosome, self.name)
            )
        self.minpos = min(self.minpos, feature.minpos)
        self.maxpos = max(self.maxpos, feature.maxpos)


def featureCmp(a, b):
    """
    General comparison function for sorting any features that have 'minpos'
    and 'maxpos' attributes.
    """
    if a.minpos == b.minpos:
        return a.maxpos - b.maxpos
    else:
        return a.minpos - b.minpos


def featureSortKey(feature):
    """Sort features by genomic interval."""
    return (feature.minpos, feature.maxpos)


def geneSortKey(gene):
    """Sort genes by interval, then id for deterministic ties."""
    return (gene.minpos, gene.maxpos, gene.id)


def featureOverlaps(a, b):
    """
    General function for determining whether feature 'a' and feature 'b' overlap.
    """
    if not (a and b):
        return False
    return intervals_overlap(a, b)


def featureContains(a, b):
    """
    General function for determining whether feature 'a' contains
    feature 'b'.  Note that both features must have 'minpos'
    and 'maxpos' attributes.
    """
    if not (a and b):
        return False
    return interval_contains(a, b)


def feature_search(features, query, lo=0, hi=None, overlap_window=8):
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


def featureSearch(features, query, lo=0, hi=None):
    """Compatibility wrapper for ``feature_search``."""
    return feature_search(features, query, lo=lo, hi=hi)


class BaseFeature(object):
    def __init__(self, featureType, start, end, chromosome, strand, attr=None):
        self.chromosome = chromosome
        self.strand = strand
        self.parent = None
        self.minpos = min(start, end)
        self.maxpos = max(start, end)
        self.featureType = featureType
        self.attributes = {} if attr is None else attr

    def acceptor(self):
        """
        Returns the location where the acceptor dimer begins: 2nt
        upstream of the exon start position.  This is 2 before the
        start on the + strand and the exact start on the - strand.
          + Example: AG|CGTATTC
          - Example: GAATACG|CT (reverse-complement)
        """
        return self.start() - 2 if self.strand == "+" else self.start()

    def __lt__(self, other: object) -> bool:
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
        )

    def contains(self, pos, strand):
        return strand == self.strand and self.minpos <= pos <= self.maxpos

    def detailString(self):
        return (
            f"Feature: {self.featureType}\nChromosome: {self.chromosome}\n"
            f"Start: {self.start()}; End: {self.end()}; Strand: '{self.strand}'"
        )

    def donor(self):
        """
        Returns the location where the donor dimer begins: at the
        end of the exon.  This is 2 nt before the start on the - strand.
          + strand example: CGTATTC|GT
          - strand example: AC|GAATACG
        """
        return self.end() if self.strand == "+" else self.end() - 2

    def end(self):
        return self.maxpos if self.strand == "+" else self.minpos

    def __eq__(self, other):
        return (
            self.minpos == other.minpos
            and self.maxpos == other.maxpos
            and self.strand == other.strand
            and self.chromosome == other.chromosome
        )

    def gffString(self, altAttributes=None):
        """Returns a GFF-formatted representation of the feature."""
        attrs = dict(self.attributes)
        if altAttributes is not None:
            attrs.update(altAttributes)

        return "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}".format(
            self.chromosome.capitalize(),
            GFF_ID,
            self.featureType,
            self.minpos,
            self.maxpos,
            self.strand,
            dictToGFF(attrs),
        )

    def gtfString(self, transcript, genePtr, exonId):
        """Returns a GTF-formatted representation of the feature."""
        attrs = dict([(k.lower(), self.attributes[k]) for k in self.attributes])
        attrs[GTF_GENE_ID] = genePtr.id
        attrs[GTF_GENE_NAME] = genePtr.name
        attrs[GTF_TRANSCRIPT] = transcript
        attrs[GTF_EXON_ID] = exonId
        try:
            source = attrs[GTF_SOURCE]
        except KeyError:
            source = GFF_ID
        return "{}\t{}\t{}\t{}\t{}\t.\t{}\t.\t{}".format(
            self.chromosome.capitalize(),
            source,
            self.featureType,
            self.minpos,
            self.maxpos,
            self.strand,
            dictToGTF(attrs),
        )

    def __hash__(self):
        return self.__str__().__hash__()

    def __len__(self):
        return self.maxpos - self.minpos + 1

    def setParent(self, id):
        self.parent = id

    def start(self):
        return self.minpos if self.strand == "+" else self.maxpos

    def __str__(self):
        return f"{self.chromosome} {self.featureType}: {self.start()}-{self.end()} ({self.strand})"


class Exon(BaseFeature):
    def __init__(self, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, RecordType.EXON, start, end, chromosome, strand, attr)
        self.parents = []

    def __str__(self):
        if self.parents:
            return (
                f"{self.featureType} {self.start()}-{self.end()}({self.strand}) "
                f"(isoforms: {','.join([x.id for x in self.parents])})"
            )
        else:
            return f"{self.featureType} {self.start()}-{self.end()}({self.strand})"

    def addParent(self, isoform):
        if isoform not in self.parents:
            self.parents.append(isoform)


class Isoform(BaseFeature):
    def __init__(self, id, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, ISOFORM_TYPE, start, end, chromosome, strand, attr)
        self.id = id
        if self.id == "ENSG00000149256":
            raise AttributeError("ISOFORM USES GENE NAME {}".format(id))
        self.features = []
        self.exons = []
        self.exonMap = {}

    def __eq__(self, other):
        return self.featureType == other.featureType and self.id == other.id

    def acceptorList(self):
        """
        Returns a list of acceptor positions for the isoform.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [e.acceptor() for e in self.exons[1:]]

    def addExon(self, exon: Exon) -> bool:
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

    def addFeature(self, feature: BaseFeature) -> None:
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

    def detailString(self):
        return (
            f"Isoform {self.id}\nStart: {self.start()}; End: {self.end()}; "
            f"Strand: {self.strand}\nExons: [{self.exonString()}]\n"
        )

    def donorList(self):
        """
        Returns a list of donor positions for the isoform.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [e.donor() for e in self.exons[:-1]]

    def exonString(self):
        return ",".join([str(e) for e in self.exons])

    def getFeatureList(self, featureType):
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
        result = []
        # Always sort in ascending order by position
        exonList = sorted(self.exons, key=featureSortKey)
        for i in range(len(exonList)):
            exon = exonList[i]
            result.append(exon.gtfString(self.id, self.parent, i + 1))
        return result

    def sortedExons(self):
        """
        Sorts the exons in an isoform based on its strand and
        returns the sorted list of exon objects.
        """
        self.exons.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return self.exons

    def sortedIntrons(self):
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

    def __str__(self):
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons)} exons/cds, "
            f"range {self.minpos} to {self.maxpos}"
        )


# Just as exons are part of an isoform, so CDS elements are part of an mRNA sequence:
class CDS(Exon):
    def __init__(self, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, RecordType.CDS, start, end, chromosome, strand, attr)
        self.parents = []

    def __lt__(self, other: object) -> bool:
        """Compare CDS records, including feature type for tied intervals."""
        if not isinstance(other, BaseFeature):
            return NotImplemented
        return (self.chromosome, self.minpos, self.maxpos, self.featureType) < (
            other.chromosome,
            other.minpos,
            other.maxpos,
            other.featureType,
        )

    def __eq__(self, o):
        """Compare CDS equality by type and genomic location."""
        return (
            self.featureType == o.featureType
            and self.minpos == o.minpos
            and self.maxpos == o.maxpos
            and self.strand == o.strand
            and self.chromosome == o.chromosome
        )

    def __str__(self):
        return (
            f"{self.featureType} {self.start()}-{self.end()} "
            f"(isoforms: {','.join([x.id for x in self.parents])})"
        )


# We treat UTR records the same way as CDS records
class FP_UTR(CDS):
    def __init__(self, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, RecordType.FIVE_PRIME_UTR, start, end, chromosome, strand, attr)
        self.parents = []


class TP_UTR(CDS):
    def __init__(self, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, RecordType.THREE_PRIME_UTR, start, end, chromosome, strand, attr)
        self.parents = []


class mRNA(Isoform):
    """
    An mRNA acts like an isoform in that it is associated with a parent gene
    and contains a number of coding sequences (CDS).
    """

    def __init__(self, id, start, end, chromosome, strand, attr=None):
        BaseFeature.__init__(self, RecordType.MRNA, start, end, chromosome, strand, attr)
        self.id = id
        self.exons = []
        self.features = []
        self.cds = []
        self.cdsMap = {}
        self.start_codon = None
        self.end_codon = None

    def acceptorList(self):
        """Returns a list of acceptor positions for the mRNA."""
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [c.acceptor() for c in self.cds[1:]]

    def addCDS(self, cds: CDS) -> bool:
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

    def donorList(self):
        """Returns a list of donor positions for the mRNA."""
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return [c.donor() for c in self.cds[:-1]]

    def endCodon(self):
        """Returns the end codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.end_codon:
            self.findCodons()
        return self.end_codon

    def findCodons(self):
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

    def inferCodons(self):
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

    def getUTRs(self):
        """Returns a list of all UTR records in the mRNA object."""
        return [
            c
            for c in self.cds
            if c.featureType in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]
        ]

    def gffStrings(self):
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
            if self.parent is None:
                raise RuntimeError("mRNA parent gene is required before writing start_codon")
            return (
                f"{self.chromosome}\t{GFF_ID}\tstart_codon\t{self.start_codon[0]}"
                f"\t{self.start_codon[1]}\t.\t{self.strand}\t.\tgene_id "
                f'"{self.parent.id}"; transcript_id "{self.id}"'
            )
        return None

    def gtfStopCodon(self) -> str | None:
        if self.end_codon:
            if self.parent is None:
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

        # Always sort in ascending order by position
        cdsList = sorted(self.cds, key=featureSortKey)
        for i in range(len(cdsList)):
            c = cdsList[i]
            result.append(c.gtfString(self.id, self.parent, i + 1))

        codonString = self.gtfStopCodon() if self.strand == "+" else self.gtfStartCodon()
        if codonString:
            result.append(codonString)

        return result

    def sortedCDS(self):
        """
        Sorts the CDS in an mRNA based on its strand and returns the sorted list of CDS objects.
        """
        self.cds.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return self.cds

    def sortedExons(self, **args):
        """
        Infers an exon list from a CDS list and sorts the list based on strand.  Usually
        this means a 5' UTR record abuts a CDS record or a CDS record abuts a 3' UTR record,
        in which case an expanded exon represents both.
        """
        minintron = getAttribute("minintron", 2, **args)
        cdsList = self.sortedCDS()
        result = []
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

    def startCodon(self):
        """Returns the start codon for this splice form as a duple of (start,end)
        positions, or None if there are none found.  Positions are relative to
        the start of the chromosome (1-based)."""
        if not self.start_codon:
            self.findCodons()
        return self.start_codon

    def __str__(self):
        return (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.cds)} exons/cds, "
            f"range {self.minpos} to {self.maxpos}"
        )


class Gene(BaseFeature):
    def __init__(self, id, note, start, end, chromosome, strand, name=None, attr=None):
        BaseFeature.__init__(self, RecordType.GENE, start, end, chromosome, strand, attr)
        self.id = id
        self.name = name if name is not None else id
        self.note = note
        self.isoforms = {}
        self.mrna = {}
        self.exons = []
        self.cds = []
        self.exonMap = {}
        self.cdsMap = {}
        self.string = ""

        # Codons for all transcripts/mRNA, one entry per transcript
        self.start_codons = {}
        self.end_codons = {}

        # All other features associated with genes in an annotation file, such as:
        #    3'/5' UTRs, mRNA, miRNA, siRNA, tRNA, rRNA, ncRNA, snRNA, snoRNA
        self.features = []

    def acceptorList(self):
        """
        Returns a list of acceptors for this gene.
        """
        acceptorSet = set()
        for k in self.isoforms.keys():
            acceptorSet.update(self.isoforms[k].acceptorList())
        for k in self.mrna.keys():
            acceptorSet.update(self.mrna[k].acceptorList())
        return sorted(list(acceptorSet), reverse=(self.strand == "-"))

    def addCDS(self, newmRNA, newCDS):
        """
        Adds a CDS to the gene if it's unique.
        Returns True if the CDS was added; false otherwise.
        """
        result = False
        cdsTuple = (newCDS.featureType, newCDS.minpos, newCDS.maxpos)
        try:
            cds = self.cdsMap[cdsTuple]
        except KeyError:
            cds = newCDS
            self.cds.append(cds)
            self.cdsMap[cdsTuple] = cds
            self.minpos = min(self.minpos, cds.minpos)
            self.maxpos = max(self.maxpos, cds.maxpos)
            result = True

        mrna = self.addmRNA(newmRNA)
        mrna.addCDS(cds)
        self.start_codons[mrna.id] = mrna.startCodon()
        self.end_codons[mrna.id] = mrna.endCodon()
        return result

    def addExon(self, newIsoform, newExon):
        """
        Adds an exon to the gene if it's unique.
        Returns True if the exon was added; false otherwise.
        """
        if newIsoform is None:
            raise ValueError("Illegal null isoform in exon {}".format(newExon))
        result = False
        posTuple = (newExon.minpos, newExon.maxpos)
        try:
            exon = self.exonMap[posTuple]
        except KeyError:
            exon = newExon
            self.exons.append(exon)
            self.exonMap[posTuple] = exon
            self.minpos = min(self.minpos, exon.minpos)
            self.maxpos = max(self.maxpos, exon.maxpos)
            result = True

        isoform = self.addIsoform(newIsoform)
        isoform.addExon(exon)
        return result

    def addFeature(self, feature: BaseFeature) -> None:
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

    def addIsoform(self, isoform):
        isoform.setParent(self)
        return self.isoforms.setdefault(isoform.id, isoform)

    def addmRNA(self, mrna):
        mrna.parent = self
        return self.mrna.setdefault(mrna.id, mrna)

    def detailString(self):
        result = (
            f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
            f"(len={len(self)}, strand={self.strand}), {len(self.exons) + len(self.cds)} "
            f"exons/cds, range {self.minpos} to {self.maxpos}\n"
        )
        result += "\n  ".join([x.detailString() for x in self.isoforms.values()])
        return result

    def donorList(self):
        """
        Returns a list of donors for this gene.
        """
        donorSet = set()
        for k in self.isoforms.keys():
            donorSet.update(self.isoforms[k].donorList())
        for k in self.mrna.keys():
            donorSet.update(self.mrna[k].donorList())
        return sorted(list(donorSet), reverse=(self.strand == "-"))

    def endCodons(self):
        """Returns a dictionary of splice forms and their associated end codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.end_codons

    def getFeatureList(self, featureType):
        return [f for f in self.features if f.featureType == featureType]

    def getIntrons(self):
        """
        Returns a list of duples containing start/end positions of introns in this gene.
        """
        result = {}
        for iid in self.isoforms.keys():
            exons = self.isoforms[iid].sortedExons()
            for i in range(1, len(exons)):
                key = (exons[i - 1].end(), exons[i].start())
                result[key] = 1

        for mid in self.mrna.keys():
            cds = self.mrna[mid].sortedExons()
            for i in range(1, len(cds)):
                key = (cds[i - 1].end(), cds[i].start())
                result[key] = 1

        return result.keys()

    def getIsoform(self, id):
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

    def gffStrings(self):
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
                allExons.sort(key=featureSortKey)

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

    def isSingleExon(self):
        return len(self.exons) == 1

    def sortedExons(self):
        """Returns a list of all exons inferred by the gene model, sorted 5' to 3'."""
        tmpset = set()
        for iid in self.isoforms.keys():
            tmpset.update(self.isoforms[iid].sortedExons())

        # Avoid returning duplicates:
        stored = set([(e.minpos, e.maxpos) for e in tmpset])
        for mid in self.mrna.keys():
            tmpset.update(
                [e for e in self.mrna[mid].sortedExons() if (e.minpos, e.maxpos) not in stored]
            )

        result = list(tmpset)
        result.sort(key=featureSortKey, reverse=(self.strand == "-"))
        return result

    def startCodons(self):
        """Returns a dictionary of splice forms and their associated start codon locations.
        A codon location is given as a duple of (start,end) positions."""
        return self.start_codons

    def __str__(self):
        if not self.string:
            self.string = (
                f"{self.id} ({self.chromosome}): {self.start()}-{self.end()} "
                f"(len={len(self)}, strand={self.strand}), {len(self.cds) + len(self.exons)} "
                f"exons/cds, range {self.minpos} to {self.maxpos}"
            )
        return self.string


class PseudoGene(Gene):
    def __init__(self, id, note, start, end, chromosome, strand, name=None, attr=None):
        BaseFeature.__init__(self, RecordType.PSEUDOGENE, start, end, chromosome, strand, attr)
        self.id = id
        self.name = name
        self.note = note
        self.features = []
        self.exons = []
        self.cds = []
        self.mrna = {}
        self.isoforms = {}
        self.exonMap = {}
        self.cdsMap = {}
        self.string = ""

        self.start_codons = {}
        self.end_codons = {}

    def detailString(self):
        return self.__str__()


class GeneModel(object):
    def __init__(self, gffPath, **args):
        """Instantiates a GeneModel object.  If a path is provided, this will
        load gene models from the given file."""
        # 2-dimensional model indexed by chromosome then gene
        self.allGenes = {}
        self.allChr = {}
        self.foundTypes = {}
        self.model = {}
        self.mRNAforms = {}
        self.sorted = {}

        if gffPath:  # Load gene models from a file
            if isinstance(gffPath, str) and not os.path.exists(gffPath):
                raise ValueError("Gene model file not found: {}".format(gffPath))

            self.loadGeneModel(gffPath, **args)

            if not self.model:
                raise ValueError("No gene models found in {}".format(gffPath))
            self.makeSortedModel()

    def __contains__(self, gene):
        """Returns true if a gene is in the model; false otherwise."""
        return str(gene) in self.allGenes

    def addChromosome(self, start, end, name):
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

    def addGene(self, gene):
        """Adds a gene to a gene model.  Raises a ValueError if the
        gene has already been added."""
        if gene.id in self.model[gene.chromosome] or str(gene) in self.allGenes:
            raise ValueError("Gene {} already stored in gene model".format(gene.id))

        self.model[gene.chromosome][gene.id] = gene
        self.allGenes[str(gene)] = gene

    def binarySearchGenes(self, geneList, loc, lo, hi, pfx=""):
        """Quasi-binary search through a gene list.  Performs binary search to a point,
        then performs linear search within a refined gene range.  Real binary search may
        not be possible since gene models may overlap or contain one within another."""
        if (hi - lo) <= GENE_SEARCH_RANGE:
            loBound = max(0, lo - GENE_SEARCH_MARGIN)
            hiBound = min(hi + GENE_SEARCH_MARGIN, len(geneList))
            for gene in geneList[loBound:hiBound]:
                if gene.contains(loc, gene.strand):
                    return (gene, gene)
            return (None, None)
        else:
            midpt = (hi + lo) // 2
            midGene = geneList[midpt]
            if midGene.contains(loc, midGene.strand):
                return (midGene, midGene)
            elif loc > midGene.maxpos:
                return self.binarySearchGenes(geneList, loc, midpt + 1, hi, pfx + " ")
            elif loc < midGene.minpos:
                return self.binarySearchGenes(geneList, loc, lo, midpt - 1, pfx + " ")

    def cleanName(self, s):
        """
        Some feature names include URL characters that we may wish to fix.
        """
        revised = unquote(s)
        revised = revised.replace(",", "")
        return revised.replace(" ", "-")

    def getAllAcceptors(self, geneFilter=defaultGeneFilter):
        """
        Returns a dictionary of all known acceptor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model.keys():
            result[chrom] = self.getKnownAcceptors(chrom, geneFilter)
        return result

    def getAllDonors(self, geneFilter=defaultGeneFilter):
        """
        Returns a dictionary of all known donor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model.keys():
            result[chrom] = self.getKnownDonors(chrom, geneFilter)
        return result

    def getAllGeneIds(self, geneFilter=defaultGeneFilter):
        """Returns a list of ids for all genes stored."""
        return [g.id for g in self.allGenes.values() if geneFilter(g)]

    def getAllGenes(self, geneFilter=defaultGeneFilter, **args):
        """Returns a list of all genes stored."""
        verbose = getAttribute("verbose", False, **args)
        indicator = ProgressIndicator(10000, verbose=verbose)
        result = []
        for g in self.allGenes.values():
            indicator.update()
            if geneFilter(g):
                result.append(g)
        indicator.finish()
        return result

    def getAnnotation(self, key, annotDict, default=None):
        """
        Convenience method for retrieving a value from an annotation dictionary
        """
        try:
            return annotDict[key]
        except KeyError:
            return default

    def getAnnotationDict(self, s):
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

    def getChromosome(self, chrName: str) -> Chromosome | None:
        """Returns a simple record with basic chromosome information."""
        try:
            return self.allChr[chrName]
        except KeyError:
            return None

    def getChromosomes(self):
        """Returns a list of all chromosomes represented in the model."""
        return self.model.keys()

    def getFeatureList(self, featureType):
        """
        Returns a list of all features of the given type found in all genes.
        """
        result = []
        for gene in self.allGenes.values():
            fList = gene.getFeatureList(featureType)
            if fList:
                result += fList
        return result

    def getGene(self, chrom: str, geneId: str) -> Gene | None:
        """
        Returns a gene from within a chromosome.  The gene will contain
        information on all exons within it.
        """
        try:
            return self.model[chrom.lower()][geneId]
        except KeyError:
            return None

    def getGeneByName(self, id: str) -> Gene | None:
        """
        Returns a gene with the given id if it exists.
        """
        for k in self.model.keys():
            try:
                return self.model[k][id.upper()]
            except KeyError:
                pass
        return None

    def getGeneFromLocations(
        self, chrom: str, startPos: int, endPos: int, strand: str
    ) -> Gene | None:
        """
        Finds the gene within the given chromosome that contains the given start
        and end positions.  Uses quasi-binary search through the sorted list of genes.
        """
        # Only get genes on the correct strand:
        try:
            geneList = self.sorted[chrom.lower()][strand]
        except KeyError:
            raise KeyError(f"Key {chrom.lower()} not found in {','.join(self.sorted.keys())}")

        (logene, higene) = self.binarySearchGenes(geneList, startPos, 0, len(geneList) - 1)

        if not (logene or higene):
            return None
        elif logene.contains(startPos, strand) or logene.contains(endPos, strand):
            return logene
        elif higene.contains(startPos, strand) or higene.contains(endPos, strand):
            return higene

        return None

    def getGeneRecords(self, chrom, geneFilter=defaultGeneFilter, **args):
        """
        Returns a list of all gene instances represented within a given chromosome.
        The gene list may be filtered by changing the geneFilter function.
        """
        verbose = getAttribute("verbose", False, **args)
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

    def getGenes(self, chrom):
        """
        Returns a list of all genes represented within a given chromosome.
        """
        try:
            return self.model[chrom.lower()].keys()
        except KeyError:
            return []

    def getGenesInRange(self, chrom, minpos, maxpos, strand=None):
        """
        Returns a list of all gene instances represented within a given chromosome
        that overlap the given range.  If no strand is specified, this will
        return all genes on both strands that overlap the range.
        """
        result = []
        for g in self.getGeneRecords(chrom):
            if g.maxpos < minpos or g.minpos > maxpos:
                continue
            if strand and g.strand != strand:
                continue
            result.append(g)
        return result

    def getKnownAcceptors(self, chrom, geneFilter=defaultGeneFilter):
        """
        Returns a dictionary of all known acceptor sites for the chromosome,
        indexed by strand.
        """
        result = {"-": set([]), "+": set([])}
        for g in self.getGeneRecords(chrom, geneFilter):
            result[g.strand].update(g.acceptorList())
        return result

    def getKnownDonors(self, chrom, geneFilter=defaultGeneFilter):
        """
        Returns a dictionary of all known donor sites for the chromosome,
        indexed by strand.
        """
        result = {"-": set([]), "+": set([])}
        for g in self.getGeneRecords(chrom, geneFilter):
            result[g.strand].update(g.donorList())
        return result

    def getParent(
        self, s: str, chrom: str, searchGenes: bool = True, searchmRNA: bool = True
    ) -> Gene | mRNA | None:
        """
        Parent identifiers are not stored in a consistent manner.  We may have
        'AT1G01160', 'AT1G01160.1' or '12345.AT1G01160' or possibly something else.
        This method looks for the most specific candidate name to identify a record's parent.
        """

        def getSubnames(fullString: str, d: str) -> list[str]:
            parts = fullString.split(d)
            result = []
            # most-to-least specific:
            # 'a-b-c-d' -> ['a-b-c', 'a-b', 'a']
            # 'a,b,c' -> ['a,b', 'a']
            for i in range(len(parts) - 1, 0, -1):
                result.append(d.join(parts[:i]))
            return result

        parentString = s.upper()
        # First try the easiest method:
        if searchmRNA:
            try:
                return self.mRNAforms[chrom][parentString]
            except KeyError:
                pass

        if searchGenes:
            try:
                return self.model[chrom][parentString]
            except KeyError:
                pass

        # Iterate over a set of known delimiters
        delim = [c for c in FORM_DELIMITERS if c in parentString]
        candidates = [parentString]
        for c in delim:
            candidates += getSubnames(parentString, c)

        # First try known mRNA records
        if searchmRNA and chrom in self.mRNAforms:
            commasep = parentString.split(",")
            additional = list(commasep)
            for c in commasep:
                additional += c.split(".")
            candidates += additional
            for c in candidates:
                try:
                    return self.mRNAforms[chrom][c]
                except KeyError:
                    pass

        # Next try known genes;
        if searchGenes and chrom in self.model:
            for c in candidates:
                if c in self.model[chrom]:
                    return self.model[chrom][c]
        return None

    def getRecordTypes(self):
        """Returns a list of all record types found in the input file."""
        return [k for k in self.foundTypes if self.foundTypes[k]]

    def isoformDict(self, **args):
        """Returns a dictionary that maps gene names to their corresponding
        isoform identifiers.  Each gene is associated with a set of isoform ids."""
        result = {}
        for g in self.getAllGenes(**args):
            result[g.id] = set(g.mrna.keys() + g.isoforms.keys())
        return result

    def loadGeneModel(self, gffRecords, **args):
        """
        Reads a tab-delimited gene annotation GFF file and stores information
        on chromosomes, the genes within each chromosome and exons within each gene.

        Parameters:
          'gffRecords'   - source of GFF records; may be a file path, a file stream
                           or a list/set of strings
          'requireNotes' - require annotations for all gene records (default=False)
          'chromosomes'  - chromosome name or list of chromosomes to store (default=all)
          'verbose'      - provide verbose feedback (default=False)
          'ignoreErrors' - ignore error conditions (default=False)
        """
        verbose = getAttribute("verbose", False, **args)
        requireNotes = getAttribute("requireNodes", False, **args)
        ignoreErrors = getAttribute("ignoreErrors", False, **args)
        chromosomes = None

        # Convenience method for handling exceptions that could be ignored:
        def conditionalException(s: str) -> None:
            if not ignoreErrors:
                raise RuntimeError(s)

        if isinstance(gffRecords, str):
            if verbose:
                sys.stderr.write("Loading and validating gene models in {}\n".format(gffRecords))
            instream = ez_open(gffRecords)
        elif isinstance(gffRecords, (list, set, tuple)):
            if verbose:
                sys.stderr.write(
                    f"Loading and validating gene models in {len(gffRecords)} records\n"
                )
            instream = gffRecords
        else:
            raise ValueError(
                "Unrecognized GFF record source "
                f"({type(gffRecords).__name__}); must be file path or a list/set of strings."
            )

        if "chromosomes" in args:
            chrParam = args["chromosomes"]
            if isinstance(chrParam, (list, tuple, set)):
                chromosomes = [s.lower() for s in chrParam]
            elif chrParam:
                chromosomes = [chrParam]
            if verbose:
                sys.stderr.write("GeneModel loading chromosomes {}\n".format(",".join(chromosomes)))

        self.model = {}
        self.mRNAforms = {}
        self.allGenes = {}
        self.foundTypes = {}
        self.allChr = {}
        lineCtr = 0

        geneAlias = {}
        geneCount = 0
        exonCount = 0
        isoCount = 0
        mrnaCount = 0
        cdsCount = 0
        badLines = 0
        indicator = ProgressIndicator(1000000, verbose=verbose)

        for line in instream:
            lineCtr += 1
            indicator.update()
            s = line.rstrip()
            if not s or s[0] == "#":
                continue  # GFF comments:

            parts = s.split("\t")

            if len(parts) < 7:
                badLines += 1
                if verbose:
                    sys.stderr.write(
                        f"line {lineCtr}: invalid GFF format "
                        "(not enough columns); file may be corrupt\n"
                    )
                if badLines >= MAX_BAD_LINES:
                    sys.stderr.write("\nInput GFF file appears to be invalid; aborting.\n")
                    raise ValueError("Invalid GFF input file")
                continue

            annots = self.getAnnotationDict(parts[-1])
            # Columns in GFF file
            # Since chromosome name and record type may be used in
            # comparisions, make them all lowercase.
            chrName = parts[0].lower()
            if chromosomes and chrName not in chromosomes:
                continue

            # Convert input record type into a known enum-backed domain and map
            # aliases (e.g. predicted_gene -> gene) before processing.
            try:
                recType = coerce_enum(parts[2].lower(), RecordType, field="record_type").value
            except ValueError as exc:
                raise ValueError(f"line {lineCtr}: unknown record type '{parts[2]}'") from exc
            recType = RECTYPE_MAP.get(recType, recType)

            startPos = int(parts[3])
            endPos = int(parts[4])
            try:
                strand = coerce_enum(parts[6], Strand, field="strand").value
            except ValueError as exc:
                raise ValueError(f"line {lineCtr}: unknown strand '{parts[6]}'") from exc

            # Track all known record types and which were stored
            self.foundTypes[recType] = recType in KNOWN_RECTYPES and recType not in IGNORE_RECTYPES
            if not self.foundTypes[recType]:
                continue

            # Many gene types possible, usually ending in '_gene'
            if recType in [RecordType.GENE, RecordType.PSEUDOGENE]:
                # Id is used for comparision, so use all uppercase
                try:
                    gid = annots[ID_FIELD].upper()
                except KeyError:
                    raise ValueError(f"line {lineCtr}: {recType} record has no ID field:\n{line}\n")

                name = self.getAnnotation(NAME_FIELD, annots, None)
                if name:
                    name = self.cleanName(name)

                note = self.getAnnotation(NOTE_FIELD, annots)
                if not note and requireNotes:
                    continue

                if strand not in VALID_STRANDS:
                    conditionalException(f"line {lineCtr}: {recType} record with unknown strand")

                if chrName not in self.model:
                    self.model[chrName] = {}
                    self.addChromosome(1, endPos, chrName)

                if recType == RecordType.PSEUDOGENE:
                    gene = PseudoGene(gid, note, startPos, endPos, chrName, strand, name, annots)
                else:  # genes and predicted_genes
                    gene = Gene(gid, note, startPos, endPos, chrName, strand, name, annots)

                try:
                    other = self.allGenes[str(gene)]
                    conditionalException(
                        f"line {lineCtr}: gene {gene.id} associated with multiple loci: "
                        f"{other.minpos}-{other.maxpos} and {startPos}-{endPos}"
                    )
                except KeyError:
                    pass

                self.allChr[chrName].update(gene)
                self.addGene(gene)
                # Map both relationships
                geneAlias[gene.name.upper()] = gene.id.upper()
                geneAlias[gene.id.upper()] = gene.name.upper()
                geneCount += 1

            elif recType in [RecordType.EXON, RecordType.PSEUDOGENIC_EXON]:
                if chrName not in self.model:
                    continue

                # Parent identifiers for exons
                parent = None
                tried = set()
                for k in POSSIBLE_GENE_FIELDS:
                    if (k not in annots) or (annots[k] in tried):
                        continue
                    isoName = annots[k]
                    parent = self.getParent(isoName, chrName)
                    if parent:
                        break
                    tried.add(isoName)
                if not parent:
                    continue

                if parent.featureType == RecordType.MRNA:
                    gene = parent.parent
                    if gene is None:
                        conditionalException(
                            f"line {lineCtr}: mRNA parent is missing gene for exon record"
                        )
                        continue
                else:
                    gene = parent

                # Next get the isoform if it already exists:
                isoform = None
                isoName = ""
                tried = set()
                # Look through known transcript ID keys for a suitable attribute:
                for k in POSSIBLE_FORM_FIELDS:
                    if k not in annots:
                        continue
                    if annots[k] in tried:
                        continue
                    isoName = annots[k]
                    if isoName in gene.isoforms:
                        isoform = gene.isoforms[isoName]
                        break
                    else:
                        tried.add(isoName)

                # If there is no form or name to trace we're out of luck:
                if not (isoform or isoName):
                    continue

                if strand in VALID_STRANDS and strand != gene.strand:
                    conditionalException(
                        f"line {lineCtr}: exon strand ({strand}) != gene "
                        f"strand ({gene.strand}) for {gene.id}"
                    )
                else:
                    strand = gene.strand

                # Must add isoform to gene before adding exon
                if not isoform:
                    isoAttr = {
                        PARENT_FIELD: gene.id,
                        NAME_FIELD: isoName,
                        ID_FIELD: isoName,
                    }
                    isoform = Isoform(isoName, startPos, endPos, chrName, strand, attr=isoAttr)
                    isoCount += 1

                exon = Exon(startPos, endPos, chrName, strand, annots)
                exonCount = exonCount + 1 if gene.addExon(isoform, exon) else exonCount

            elif recType in [RecordType.MRNA, RecordType.PSEUDOGENIC_TRANSCRIPT]:
                if chrName not in self.model:
                    conditionalException(
                        f"line {lineCtr}: mRNA with missing chromosome dictionary {chrName} "
                        f"(known: {','.join(self.model.keys())})"
                    )

                if ID_FIELD not in annots:
                    conditionalException(f"line {lineCtr}: mRNA with missing ID")

                id = annots[ID_FIELD].upper()
                parent = self.getAnnotation(PARENT_FIELD, annots)
                if not parent:
                    continue

                parent = parent.upper()
                gene = self.getParent(parent, chrName)
                if not gene:
                    try:
                        alias = geneAlias[parent]
                        gene = self.getParent(alias, chrName)
                    except KeyError:
                        if verbose:
                            if alias == parent:
                                sys.stderr.write(
                                    f"line {lineCtr}: no gene {parent} found for {recType}\n"
                                )
                            else:
                                sys.stderr.write(
                                    f"line {lineCtr}: no gene '{parent}' or "
                                    f"'{alias}' found for {recType}\n"
                                )
                        continue

                if strand in VALID_STRANDS and strand != gene.strand:
                    conditionalException(
                        f"line {lineCtr}: mRNA strand ({strand}) does not "
                        f"match gene strand ({gene.strand})"
                    )
                else:
                    strand = gene.strand

                mrnaAttr = {PARENT_FIELD: gene.id, NAME_FIELD: id, ID_FIELD: id}
                mrna = mRNA(id, startPos, endPos, chrName, strand, attr=mrnaAttr)
                gene.addmRNA(mrna)
                self.mRNAforms[chrName] = self.mRNAforms.setdefault(chrName, {})
                self.mRNAforms[chrName][id] = mrna
                mrnaCount += 1

            elif recType in CDS_TYPES:
                if chrName not in self.model:
                    conditionalException(
                        f"line {lineCtr}: {recType} has unrecognized chromosome: {chrName} "
                        f"(known: {','.join(self.model.keys())})"
                    )
                if chrName not in self.mRNAforms:
                    conditionalException(
                        f"line {lineCtr}: {recType} has unrecognized chromosome: {chrName} "
                        f"(known: {','.join(self.mRNAforms.keys())})"
                    )

                mrna = self.getParent(annots[PARENT_FIELD], chrName, searchGenes=False)
                if not mrna:
                    if verbose:
                        sys.stderr.write(
                            f"line {lineCtr}: no mRNA {annots[PARENT_FIELD]} found for {recType}\n"
                        )
                    continue
                if mrna.featureType != RecordType.MRNA:
                    conditionalException(
                        f"line {lineCtr}: parent {annots[PARENT_FIELD]} is not an mRNA record"
                    )
                    continue
                gene = mrna.parent
                if gene is None:
                    conditionalException(
                        f"line {lineCtr}: mRNA {mrna.id} is missing a parent gene for CDS record"
                    )
                    continue

                if strand in VALID_STRANDS and strand != mrna.strand:
                    conditionalException(
                        f"line {lineCtr}: CDS strand ({strand}) does not "
                        f"match mRNA strand ({mrna.strand})"
                    )
                else:
                    strand = mrna.strand

                cds = cdsFactory(recType, startPos, endPos, chrName, strand, annots)

                if gene.addCDS(mrna, cds):
                    cdsCount += 1

            elif recType == RecordType.CHROMOSOME:
                self.addChromosome(startPos, endPos, chrName)

            elif recType in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]:
                if chrName not in self.model:
                    continue
                if PARENT_FIELD not in annots:
                    continue

                parentId = annots[PARENT_FIELD]
                parent = self.getParent(parentId, chrName)
                if parent:
                    if strand in VALID_STRANDS and strand != parent.strand:
                        conditionalException(
                            f"line {lineCtr}: {recType} strand ({strand}) does not match "
                            f"parent strand ({parent.strand})"
                        )
                    else:
                        strand = parent.strand
                    parent.addFeature(
                        BaseFeature(recType, startPos, endPos, chrName, strand, annots)
                    )
                else:
                    if verbose:
                        sys.stderr.write(
                            f"line {lineCtr}: no parent {parentId} found for {recType}\n"
                        )
                    continue

            elif recType not in IGNORE_RECTYPES:
                if chrName not in self.model:
                    continue
                if PARENT_FIELD not in annots:
                    continue

                parent = annots[PARENT_FIELD].split(".")[0]
                if parent not in self.model[chrName]:
                    continue

                gene = self.model[chrName][parent]
                try:
                    gene.addFeature(BaseFeature(recType, startPos, endPos, chrName, strand, annots))
                except ValueError as e:
                    conditionalException(f"line {lineCtr}: {e}")

        indicator.finish()
        if verbose:
            if geneCount > 0:
                sys.stderr.write(
                    "Loaded {} genes with {} isoforms, {} exons (avg. {:.1f}/gene), "
                    "{} mRNA, {} CDS (avg. {:.1f}/gene)\n".format(
                        comma_format(geneCount),
                        comma_format(isoCount),
                        comma_format(exonCount),
                        float(exonCount) / geneCount,
                        comma_format(mrnaCount),
                        comma_format(cdsCount),
                        float(cdsCount / geneCount),
                    )
                )
            else:
                sys.stderr.write("** Warning: no genes loaded!\n")

    def makeSortedModel(self):
        self.sorted = {}
        for chrom in self.model.keys():
            self.sorted[chrom] = {}
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

    def writeGFF(self, gffPath, **args):
        """Writes a complete gene model out to a GFF file."""
        geneFilter = getAttribute("geneFilter", defaultGeneFilter, **args)
        geneSubset = getAttribute("geneSet", None, **args)
        verbose = getAttribute("verbose", False, **args)

        outStream = open(gffPath, "w") if isinstance(gffPath, str) else gffPath
        chromList = sorted(self.allChr.keys())
        indicator = ProgressIndicator(10000, verbose=verbose)
        for c in chromList:
            chrom = self.getChromosome(c)
            if chrom is None:
                continue
            outStream.write("{}\n".format(chrom.gffString()))
            genes = self.getGeneRecords(c, geneFilter)
            if geneSubset:
                genes = [g for g in genes if g.id in geneSubset or g.name in geneSubset]
            genes.sort(key=geneSortKey)
            for g in genes:
                indicator.update()
                strings = g.gffStrings()
                if strings:
                    outStream.write("{}\n".format(strings))
        indicator.finish()

    def writeGTF(self, gtfPath, geneFilter=defaultGeneFilter, **args):
        """Writes a complete gene model out to a GTF file."""
        verbose = getAttribute("verbose", False, **args)
        outStream = open(gtfPath, "w") if isinstance(gtfPath, str) else gtfPath
        chromList = sorted(self.allChr.keys())
        indicator = ProgressIndicator(10000, verbose=verbose)
        for c in chromList:
            genes = self.getGeneRecords(c, geneFilter)
            genes.sort(key=geneSortKey)
            for g in genes:
                indicator.update()
                strings = g.gtfStrings()
                if strings:
                    outStream.write("{}\n".format(strings))
        indicator.finish()
