"""Canonical enum domains used across SpliceGrapher."""

from __future__ import annotations

from enum import StrEnum


class Strand(StrEnum):
    PLUS = "+"
    MINUS = "-"
    UNKNOWN = "."


class RecordType(StrEnum):
    GENE = "gene"
    PSEUDOGENE = "pseudogene"
    PREDICTED_GENE = "predicted_gene"
    MRNA = "mrna"
    PSEUDOGENIC_TRANSCRIPT = "pseudogenic_transcript"
    EXON = "exon"
    PSEUDOGENIC_EXON = "pseudogenic_exon"
    CDS = "cds"
    CDS_PREDICTED = "cds_predicted"
    FIVE_PRIME_UTR = "five_prime_utr"
    THREE_PRIME_UTR = "three_prime_utr"
    INTRON = "intron"
    CHROMOSOME = "chromosome"
    PROTEIN = "protein"
    NONUNIQUE = "nonunique"
    MRNA_TE_GENE = "mRNA_TE_gene"
    TRANS_ELE_GENE = "transposable_element_gene"
    GRAPH = "graph"
    CLUSTER = "cluster"
    PARENT = "parent"
    CHILD = "child"


class NodeDisposition(StrEnum):
    KNOWN = "known"
    PREDICTED = "predicted"
    UNRESOLVED = "unresolved"


class AlternativeSplicingEvent(StrEnum):
    IR = "IR"
    ALT5 = "A5"
    ALT3 = "A3"
    ALTB3 = "AB3"
    ALTB5 = "AB5"
    ES = "SE"
    ALTI = "AI"
    ALTT = "AT"


class AlternativeSplicingEventName(StrEnum):
    IR = "Intron Retention"
    ALT5 = "Alt. 5'"
    ALT3 = "Alt. 3'"
    ALTB3 = "Alt. B3"
    ALTB5 = "Alt. B5"
    ES = "Skipped Exon"
    ALTI = "Alt. Init."
    ALTT = "Alt. Term."


ALT_SPLICE_EVENT_NAME_BY_CODE: dict[
    AlternativeSplicingEvent,
    AlternativeSplicingEventName,
] = {
    AlternativeSplicingEvent.IR: AlternativeSplicingEventName.IR,
    AlternativeSplicingEvent.ALT5: AlternativeSplicingEventName.ALT5,
    AlternativeSplicingEvent.ALT3: AlternativeSplicingEventName.ALT3,
    AlternativeSplicingEvent.ALTB3: AlternativeSplicingEventName.ALTB3,
    AlternativeSplicingEvent.ALTB5: AlternativeSplicingEventName.ALTB5,
    AlternativeSplicingEvent.ES: AlternativeSplicingEventName.ES,
    AlternativeSplicingEvent.ALTI: AlternativeSplicingEventName.ALTI,
    AlternativeSplicingEvent.ALTT: AlternativeSplicingEventName.ALTT,
}

ALT_SPLICE_EVENT_CODE_BY_NAME: dict[
    AlternativeSplicingEventName,
    AlternativeSplicingEvent,
] = {name: code for code, name in ALT_SPLICE_EVENT_NAME_BY_CODE.items()}


class EdgeType(StrEnum):
    PARENT = "parent"
    CHILD = "child"


class AttrKey(StrEnum):
    ID = "ID"
    NAME = "Name"
    PARENT = "Parent"
    NOTE = "Note"
    ALT_FORM = "AltForm"
    START_CODON = "StartCodon"
    END_CODON = "EndCodon"
    ISOFORMS = "Isoforms"
    DISPOSITION = "disposition"
    EXPANDED = "Expanded"


class ShortReadCode(StrEnum):
    CHROM = "C"
    DEPTH = "D"
    READ = "R"
    JUNCTION = "J"
    SPLICE = "S"


class JunctionCode(StrEnum):
    KNOWN = "K"
    UNKNOWN = "U"
    PREDICTED = "P"
    UNLABELED = ""


class SamHeaderTag(StrEnum):
    HD = "HD"
    LN = "LN"
    SN = "SN"
    SO = "SO"
    SQ = "SQ"
    VN = "VN"


class SamHeaderLine(StrEnum):
    SQ = "@SQ"
