"""Canonical enum domains used across SpliceGrapher."""

from __future__ import annotations

from enum import Enum


class Strand(str, Enum):
    PLUS = "+"
    MINUS = "-"
    UNKNOWN = "."


class RecordType(str, Enum):
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


class NodeDisposition(str, Enum):
    KNOWN = "known"
    PREDICTED = "predicted"
    UNRESOLVED = "unresolved"


class EdgeType(str, Enum):
    PARENT = "parent"
    CHILD = "child"


class AttrKey(str, Enum):
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
