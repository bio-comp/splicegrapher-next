"""Shared constants for the splice-graph core."""

from __future__ import annotations

from SpliceGrapher.core.enums import (
    ALT_SPLICE_EVENT_CODE_BY_NAME,
    AlternativeSplicingEvent,
    AlternativeSplicingEventName,
    AttrKey,
    NodeDisposition,
    RecordType,
    Strand,
)

SOURCE_NAME = "SpliceGraph"
DEPRECATED_GENE_REC = RecordType.CLUSTER.value
GENE_REC = RecordType.GRAPH.value
PARENT_REC = RecordType.PARENT.value
CHILD_REC = RecordType.CHILD.value
VALID_GENES = frozenset({DEPRECATED_GENE_REC, GENE_REC})
VALID_RECTYPES = frozenset({DEPRECATED_GENE_REC, GENE_REC, PARENT_REC, CHILD_REC})

PARENT_ATTR = AttrKey.PARENT.value
ID_ATTR = AttrKey.ID.value
KNOWN_ATTRS = frozenset({PARENT_ATTR, ID_ATTR})
AS_KEY = AttrKey.ALT_FORM.value
START_CODON_KEY = AttrKey.START_CODON.value
END_CODON_KEY = AttrKey.END_CODON.value
ISO_KEY = AttrKey.ISOFORMS.value
DISPOSITION_KEY = AttrKey.DISPOSITION.value
KNOWN_NODE = NodeDisposition.KNOWN.value
PREDICTED_NODE = NodeDisposition.PREDICTED.value
UNRESOLVED_NODE = NodeDisposition.UNRESOLVED.value
PUTATIVE_PARENTS = "putative_parents"
PUTATIVE_CHILDREN = "putative_children"
ACCEPTORS_KEY = "acceptors"
DONORS_KEY = "donors"
ALT_MODEL_KEY = AttrKey.EXPANDED.value
VALID_STRANDS = frozenset({strand.value for strand in Strand})

NodeAttributeValue = str | set[int]


def coerce_alt_splicing_event(
    form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
) -> AlternativeSplicingEvent:
    if isinstance(form, AlternativeSplicingEvent):
        return form
    if isinstance(form, AlternativeSplicingEventName):
        return ALT_SPLICE_EVENT_CODE_BY_NAME[form]

    normalized = form.strip()
    if not normalized:
        raise ValueError("Alternative splicing form cannot be empty")

    try:
        return AlternativeSplicingEvent(normalized)
    except ValueError:
        try:
            return ALT_SPLICE_EVENT_CODE_BY_NAME[AlternativeSplicingEventName(normalized)]
        except ValueError as exc:
            raise ValueError(f"Unknown alternative splicing form: {form}") from exc


__all__ = [
    "ACCEPTORS_KEY",
    "ALT_MODEL_KEY",
    "AS_KEY",
    "CHILD_REC",
    "DEPRECATED_GENE_REC",
    "DISPOSITION_KEY",
    "DONORS_KEY",
    "END_CODON_KEY",
    "GENE_REC",
    "ID_ATTR",
    "ISO_KEY",
    "KNOWN_ATTRS",
    "KNOWN_NODE",
    "NodeAttributeValue",
    "PARENT_ATTR",
    "PARENT_REC",
    "PREDICTED_NODE",
    "PUTATIVE_CHILDREN",
    "PUTATIVE_PARENTS",
    "SOURCE_NAME",
    "START_CODON_KEY",
    "UNRESOLVED_NODE",
    "VALID_GENES",
    "VALID_RECTYPES",
    "VALID_STRANDS",
    "coerce_alt_splicing_event",
]
