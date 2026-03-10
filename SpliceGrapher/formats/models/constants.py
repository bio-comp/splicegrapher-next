"""Shared gene-model constants, type aliases, and normalization helpers."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TypeVar

from SpliceGrapher.core.enums import AttrKey, RecordType, Strand

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

GTF_GENE_ID = "gene_id"
GTF_GENE_NAME = "gene_name"
GTF_TRANSCRIPT = "transcript_id"
GTF_TRANSNAME = "transcript_name"
GTF_SOURCE = "gene_biotype"
GTF_EXON_ID = "exon_number"
GTF_PROTEIN_ID = "protein_id"

RECTYPE_MAP = {k: k for k in KNOWN_RECTYPES}
RECTYPE_MAP[RecordType.PREDICTED_GENE] = RecordType.GENE
RECTYPE_MAP[RecordType.CDS_PREDICTED] = RecordType.CDS

ISOFORM_TYPE = "isoform"

ID_FIELD = AttrKey.ID
NAME_FIELD = AttrKey.NAME
NOTE_FIELD = AttrKey.NOTE
PARENT_FIELD = AttrKey.PARENT

POSSIBLE_GENE_FIELDS = [PARENT_FIELD, GTF_GENE_ID, GTF_GENE_NAME]
POSSIBLE_FORM_FIELDS = [PARENT_FIELD, GTF_TRANSCRIPT]

GFF_ID = "SpliceGrapher"
MAX_BAD_LINES = 3

VALID_STRANDS = {strand.value for strand in Strand}

FORM_DELIMITERS = [".", "-", "_", ","]
SPLICE_DIMER_OFFSET = 2
GTF_ORDER_POLICY = "genomic_ascending"

GffRecordSource = str | list[str] | set[str] | tuple[str, ...]
AttrValue = str
AttrMap = Mapping[str | AttrKey, AttrValue]
ExtraAttrMap = Mapping[str, AttrValue] | Mapping[AttrKey, AttrValue]
AttrT = TypeVar("AttrT")


def normalize_attributes(attr: AttrMap | None) -> dict[str, str]:
    """Normalize attribute keys and values to builtin strings."""

    return {} if attr is None else {str(key): str(value) for key, value in attr.items()}
