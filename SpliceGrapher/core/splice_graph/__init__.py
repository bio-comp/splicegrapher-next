"""NetworkX-backed splice graph core package."""

from .constants import (
    AS_KEY,
    CHILD_REC,
    END_CODON_KEY,
    GENE_REC,
    ID_ATTR,
    ISO_KEY,
    KNOWN_ATTRS,
    PARENT_ATTR,
    PARENT_REC,
    SOURCE_NAME,
    START_CODON_KEY,
    VALID_GENES,
    VALID_RECTYPES,
)
from .graph import SpliceGraph
from .node import SpliceGraphNode

__all__ = [
    "AS_KEY",
    "CHILD_REC",
    "END_CODON_KEY",
    "GENE_REC",
    "ID_ATTR",
    "ISO_KEY",
    "KNOWN_ATTRS",
    "PARENT_ATTR",
    "PARENT_REC",
    "SOURCE_NAME",
    "START_CODON_KEY",
    "SpliceGraph",
    "SpliceGraphNode",
    "VALID_GENES",
    "VALID_RECTYPES",
]
