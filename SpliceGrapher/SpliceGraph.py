"""NetworkX-backed splice graph core and parser."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from dataclasses import dataclass, field
from typing import TextIO

import networkx as nx
import structlog

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import (
    ALT_SPLICE_EVENT_CODE_BY_NAME,
    AlternativeSplicingEvent,
    AlternativeSplicingEventName,
    AttrKey,
    NodeDisposition,
    RecordType,
    Strand,
)
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.progress import ProgressIndicator

LOGGER = structlog.get_logger(__name__)

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


def _coerce_alt_splicing_event(
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


@dataclass(slots=True)
class NullNode:
    """Coordinate-only node used for equality lookups by interval."""

    start: int
    end: int
    minpos: int = field(init=False)
    maxpos: int = field(init=False)

    def __post_init__(self) -> None:
        self.minpos = min(self.start, self.end)
        self.maxpos = max(self.start, self.end)


@dataclass(slots=True, init=False, eq=False)
class SpliceGraphNode:
    """Pure data container for genomic coordinates and node attributes."""

    id: str
    chromosome: str
    strand: str
    minpos: int
    maxpos: int
    start: int
    end: int
    attrs: dict[str, NodeAttributeValue]
    altFormSet: set[AlternativeSplicingEvent]
    isoformSet: set[str]
    origStart: int
    origEnd: int

    def __init__(
        self,
        id: str,
        start: int,
        end: int,
        strand: str | Strand,
        chrom: str,
    ) -> None:
        self.id = id
        self.chromosome = chrom
        self.strand = coerce_enum(strand, Strand, field="strand").value
        self.minpos = min(start, end)
        self.maxpos = max(start, end)
        if self.strand == Strand.MINUS.value:
            self.start = self.maxpos
            self.end = self.minpos
        else:
            self.start = self.minpos
            self.end = self.maxpos
        self.attrs = {}
        self.altFormSet = set()
        self.isoformSet = set()
        self.origStart = self.start
        self.origEnd = self.end

    def acceptorEnd(self) -> int:
        return self.start

    def donorEnd(self) -> int:
        return self.end

    def _sync_alt_form_attr(self) -> None:
        self.attrs[AS_KEY] = ",".join(sorted(form.value for form in self.altFormSet))

    def addAltForm(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        if isinstance(form, str) and not form.strip():
            return
        event = _coerce_alt_splicing_event(form)
        self.altFormSet.add(event)
        self._sync_alt_form_attr()

    def removeAltForm(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        if isinstance(form, str) and not form.strip():
            return
        event = _coerce_alt_splicing_event(form)
        self.altFormSet.discard(event)
        self._sync_alt_form_attr()

    def addAttribute(self, key: str, value: NodeAttributeValue) -> None:
        self.attrs[key] = value

    def addCodon(self, codon: tuple[int, int], codon_type: str) -> None:
        if len(codon) != 2:
            raise ValueError(f"Codons must be 2-tuples; received {codon!r}")
        pos = min(codon) if self.strand == Strand.PLUS.value else max(codon)
        if self.contains(pos):
            existing = self.attrs.setdefault(codon_type, set())
            if not isinstance(existing, set):
                raise TypeError(f"Codon attribute {codon_type} is not a set")
            existing.add(pos)

    def addStartCodon(self, codon: tuple[int, int]) -> None:
        self.addCodon(codon, START_CODON_KEY)

    def addEndCodon(self, codon: tuple[int, int]) -> None:
        self.addCodon(codon, END_CODON_KEY)

    def addFormsFromString(self, forms: str) -> None:
        for form in forms.split(","):
            self.addAltForm(form.strip())

    def addIsoform(self, isoform: str) -> None:
        if isoform is None:
            raise ValueError(f"Received illegal isoform for {self.id}")
        self.isoformSet.add(isoform)
        self.attrs[ISO_KEY] = ",".join(sorted(self.isoformSet))

    def addIsoformString(self, isoform_string: str) -> None:
        for isoform in isoform_string.split(","):
            cleaned = isoform.strip()
            if cleaned:
                self.addIsoform(cleaned)

    def altForms(self) -> list[str]:
        return [form.value for form in sorted(self.altFormSet, key=lambda form: form.value)]

    def altFormString(self) -> str:
        value = self.attrs.get(AS_KEY)
        return value if isinstance(value, str) else ""

    def attributeString(self) -> str:
        attr_parts: list[str] = []
        for key in sorted(self.attrs):
            if key == AS_KEY and not self.altFormSet:
                continue
            if key == ISO_KEY and not self.isoformSet:
                continue
            value = self.attrs[key]
            if isinstance(value, set):
                rendered = ",".join(str(item) for item in sorted(value))
            else:
                rendered = value
            attr_parts.append(f"{key}={rendered}")
        return ";".join(attr_parts)

    def codons(self, codon_type: str) -> list[int]:
        value = self.attrs.get(codon_type)
        if isinstance(value, set):
            return sorted(value)
        return []

    def codonString(self, codon_type: str) -> str | None:
        codons = self.codons(codon_type)
        if not codons:
            return None
        return ",".join(str(value) for value in codons)

    def contains(self, pos: int) -> bool:
        return self.minpos <= pos <= self.maxpos

    def downstreamOf(self, pos: int) -> bool:
        return self.minpos > pos if self.strand == Strand.PLUS.value else self.maxpos < pos

    def endCodons(self) -> list[int]:
        return self.codons(END_CODON_KEY)

    def endCodonString(self) -> str | None:
        return self.codonString(END_CODON_KEY)

    def hasAS(self) -> bool:
        return bool(self.altFormSet)

    def hasDisposition(self, disposition: str) -> bool:
        value = self.attrs.get(DISPOSITION_KEY)
        return value == disposition

    def isAltAcceptor(self) -> bool:
        return AlternativeSplicingEvent.ALT3 in self.altFormSet

    def isAltDonor(self) -> bool:
        return AlternativeSplicingEvent.ALT5 in self.altFormSet

    def isKnown(self) -> bool:
        return self.hasDisposition(KNOWN_NODE)

    def isPredicted(self) -> bool:
        return self.hasDisposition(PREDICTED_NODE)

    def isRetainedIntron(self) -> bool:
        return AlternativeSplicingEvent.IR in self.altFormSet

    def isSkippedExon(self) -> bool:
        return AlternativeSplicingEvent.ES in self.altFormSet

    def isUnresolved(self) -> bool:
        return self.hasDisposition(UNRESOLVED_NODE)

    def isoformList(self) -> list[str]:
        return sorted(self.isoformSet)

    def isoformString(self) -> str | None:
        value = self.attrs.get(ISO_KEY)
        return value if isinstance(value, str) else None

    def putativeChildren(self) -> set[str]:
        value = self.attrs.get(PUTATIVE_CHILDREN)
        if isinstance(value, str) and value:
            return set(value.split(","))
        return set()

    def putativeParents(self) -> set[str]:
        value = self.attrs.get(PUTATIVE_PARENTS)
        if isinstance(value, str) and value:
            return set(value.split(","))
        return set()

    def startCodons(self) -> list[int]:
        return self.codons(START_CODON_KEY)

    def startCodonString(self) -> str | None:
        return self.codonString(START_CODON_KEY)

    def update(self, minpos: int, maxpos: int) -> None:
        self.minpos = minpos
        self.maxpos = maxpos
        if self.strand == Strand.MINUS.value:
            self.start = maxpos
            self.end = minpos
        else:
            self.start = minpos
            self.end = maxpos

    def upstreamOf(self, pos: int) -> bool:
        return self.maxpos < pos if self.strand == Strand.PLUS.value else self.minpos > pos

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, (NullNode, SpliceGraphNode)):
            return NotImplemented
        return self.minpos == other.minpos and self.maxpos == other.maxpos

    def __lt__(self, other: SpliceGraphNode) -> bool:
        return (self.minpos, self.maxpos, self.id) < (other.minpos, other.maxpos, other.id)

    def __hash__(self) -> int:
        return hash((self.minpos, self.maxpos))

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1

    def __repr__(self) -> str:
        return f"{self.id} {self.start}-{self.end}"

    def __str__(self) -> str:
        return f"{self.id} {self.start}-{self.end}"


class SpliceGraph:
    """NetworkX-backed manager for splice graph topology."""

    def __init__(self, name: str, chromosome: str, strand: str | Strand) -> None:
        self.chromosome = chromosome
        self.strand = coerce_enum(strand, Strand, field="strand").value
        self.attrs: dict[str, str] = {}
        self.minpos = sys.maxsize
        self.maxpos = 0
        self._nx_graph: nx.DiGraph = nx.DiGraph()
        self.setName(name)

    @property
    def nodeDict(self) -> dict[str, SpliceGraphNode]:
        return {
            node_id: node_data["data"] for node_id, node_data in self._nx_graph.nodes(data=True)
        }

    def _sorted_nodes(self, node_ids: Iterable[str]) -> list[SpliceGraphNode]:
        nodes = [self._nx_graph.nodes[node_id]["data"] for node_id in node_ids]
        return sorted(nodes, key=lambda node: (node.minpos, node.maxpos, node.id))

    def _find_existing_node(self, start: int, end: int) -> SpliceGraphNode | None:
        probe = NullNode(start, end)
        for node in self.nodeDict.values():
            if node == probe:
                return node
        return None

    def _recompute_bounds(self) -> None:
        nodes = list(self.nodeDict.values())
        if not nodes:
            self.minpos = sys.maxsize
            self.maxpos = 0
            return
        self.minpos = min(node.minpos for node in nodes)
        self.maxpos = max(node.maxpos for node in nodes)

    def _node(self, node_or_id: str | SpliceGraphNode) -> SpliceGraphNode:
        if isinstance(node_or_id, SpliceGraphNode):
            return node_or_id
        return self._nx_graph.nodes[node_or_id]["data"]

    def addNode(self, newId: str, start: int, end: int) -> SpliceGraphNode:
        existing = self._find_existing_node(start, end)
        if existing is not None:
            return existing
        node = SpliceGraphNode(newId, start, end, self.strand, self.chromosome)
        self._nx_graph.add_node(node.id, data=node)
        self.minpos = min(self.minpos, node.minpos)
        self.maxpos = max(self.maxpos, node.maxpos)
        return node

    def addCodons(self, codon_list: Iterable[tuple[int, int]], codon_type: str) -> None:
        for codon in codon_list:
            for node in self.nodeDict.values():
                node.addCodon(codon, codon_type)

    def addEdge(self, pid: str, cid: str) -> None:
        if pid not in self._nx_graph:
            raise ValueError(f"Error adding edge: parent node {pid} not found in graph")
        if cid not in self._nx_graph:
            raise ValueError(f"Error adding edge: child node {cid} not found in graph")
        self._nx_graph.add_edge(pid, cid)

    def addEndCodons(self, codon_list: Iterable[tuple[int, int]]) -> None:
        self.addCodons(codon_list, END_CODON_KEY)

    def addStartCodons(self, codon_list: Iterable[tuple[int, int]]) -> None:
        self.addCodons(codon_list, START_CODON_KEY)

    def adjust(self, adjustment: int) -> None:
        for node in self.nodeDict.values():
            node.update(node.minpos + adjustment, node.maxpos + adjustment)
        if self.nodeDict:
            self.minpos += adjustment
            self.maxpos += adjustment

    def attributeString(self) -> str:
        return ";".join(f"{key}={self.attrs[key]}" for key in sorted(self.attrs))

    def deleteNode(self, node_id: str | SpliceGraphNode) -> SpliceGraphNode:
        node = self._node(node_id)
        self._nx_graph.remove_node(node.id)
        self._recompute_bounds()
        return node

    def getLeaves(self) -> list[SpliceGraphNode]:
        return self._sorted_nodes(
            node_id for node_id, degree in self._nx_graph.out_degree() if degree == 0
        )

    def getName(self) -> str:
        return self.attrs[ID_ATTR]

    def getNode(self, start: int, end: int) -> SpliceGraphNode | None:
        return self._find_existing_node(start, end)

    def getRoots(self) -> list[SpliceGraphNode]:
        return self._sorted_nodes(
            node_id for node_id, degree in self._nx_graph.in_degree() if degree == 0
        )

    def isEmpty(self) -> bool:
        return self._nx_graph.number_of_nodes() == 0

    def predecessors(self, node_or_id: str | SpliceGraphNode) -> list[SpliceGraphNode]:
        node = self._node(node_or_id)
        return self._sorted_nodes(self._nx_graph.predecessors(node.id))

    def predecessor_ids(self, node_id: str) -> list[str]:
        return sorted(self._nx_graph.predecessors(node_id))

    def resolvedNodes(self) -> list[SpliceGraphNode]:
        return [node for node in self.nodeDict.values() if not node.isUnresolved()]

    def setName(self, name: str) -> None:
        self.attrs[ID_ATTR] = name

    def successors(self, node_or_id: str | SpliceGraphNode) -> list[SpliceGraphNode]:
        node = self._node(node_or_id)
        return self._sorted_nodes(self._nx_graph.successors(node.id))

    def successor_ids(self, node_id: str) -> list[str]:
        return sorted(self._nx_graph.successors(node_id))

    def unresolvedNodes(self) -> list[SpliceGraphNode]:
        return [node for node in self.nodeDict.values() if node.isUnresolved()]

    def validate(self, halt: bool = False) -> str | None:
        reason: str | None = None
        if self._nx_graph.number_of_nodes() == 0:
            reason = "Graph is empty."
        else:
            roots = self.getRoots()
            leaves = self.getLeaves()
            if not roots or not leaves:
                reason = (
                    f"Graph is missing roots or leaves ({len(roots)} roots, {len(leaves)} leaves)"
                )
            elif not nx.is_directed_acyclic_graph(self._nx_graph):
                reason = "Graph contains cycles (not a valid DAG)."

        if reason and halt:
            LOGGER.error("illegal_graph_validation_failed", reason=reason)
            raise ValueError(f"Illegal graph:\n{reason}")
        return reason

    def __len__(self) -> int:
        return self.maxpos - self.minpos + 1 if self.maxpos >= self.minpos else 0

    def __str__(self) -> str:
        return (
            f"{self.getName()} ({self.strand}) {self.minpos}-{self.maxpos} "
            f"({self._nx_graph.number_of_nodes()} nodes)"
        )


class SpliceGraphParser:
    """Parse splice graph GFF files into networkx-backed ``SpliceGraph`` objects."""

    def __init__(self, fileRef: str | TextIO, *, verbose: bool = False) -> None:
        self.verbose = verbose
        self.instream = ez_open(fileRef) if isinstance(fileRef, str) else fileRef
        if self.instream is None:
            raise ValueError("No input file stream given.")
        self.graphDict: dict[str, SpliceGraph] = {}
        self._graph_keys: list[str] = []
        self.loadFromFile()

    def __iter__(self) -> SpliceGraphParser:
        return self

    def __next__(self) -> SpliceGraph:
        if self.graphId >= len(self._graph_keys):
            raise StopIteration
        key = self._graph_keys[self.graphId]
        self.graphId += 1
        return self.graphDict[key]

    def __len__(self) -> int:
        return len(self.graphDict)

    def loadFromFile(self) -> None:
        line_no = 0
        graph: SpliceGraph | None = None
        alias: dict[str, str] = {}
        edges: set[tuple[str, str]] = set()
        indicator = ProgressIndicator(100000, verbose=self.verbose)
        try:
            for line in self.instream:
                indicator.update()
                line_no += 1
                if line.startswith("#"):
                    continue
                record = line.strip()
                parts = record.split("\t")
                try:
                    rec_type = coerce_enum(parts[2].lower(), RecordType, field="record_type").value
                    start = int(parts[3])
                    end = int(parts[4])
                except (IndexError, ValueError) as exc:
                    raise ValueError(
                        f"Illegal record in splice graph file at line {line_no}:\n\t{record}"
                    ) from exc

                if rec_type not in VALID_RECTYPES:
                    raise ValueError(
                        f"Illegal record type in splice graph file at line {line_no}:\n\t{record}"
                    )

                attrs = self._parse_attributes(parts[-1], line_no)
                node_id = attrs.get(ID_ATTR)
                if node_id is None:
                    raise ValueError(
                        f"GFF attribute field '{parts[-1]}' has no ID at line {line_no}"
                    )

                if rec_type in VALID_GENES:
                    if graph is not None:
                        for parent_id, child_id in edges:
                            graph.addEdge(alias[parent_id], alias[child_id])
                    graph = SpliceGraph(name=node_id, chromosome=parts[0], strand=parts[6])
                    graph.minpos = min(start, end)
                    graph.maxpos = max(start, end)
                    for key, value in attrs.items():
                        if key not in KNOWN_ATTRS:
                            graph.attrs[key] = value
                    self.graphDict[node_id] = graph
                    alias = {}
                    edges = set()
                    continue

                if graph is None:
                    raise ValueError(f"Graph feature found before graph header at line {line_no}")

                node = graph.addNode(node_id, start, end)
                alias[node_id] = node.id
                for key, value in attrs.items():
                    if key == AS_KEY:
                        node.addFormsFromString(value)
                    elif key == ISO_KEY and value:
                        node.addIsoformString(value)
                    elif key in {START_CODON_KEY, END_CODON_KEY} and value:
                        node.attrs[key] = {int(item) for item in value.split(",")}
                    elif key not in KNOWN_ATTRS:
                        node.addAttribute(key, value)

                if PARENT_ATTR in attrs:
                    for parent in attrs[PARENT_ATTR].split(","):
                        edges.add((parent, node_id))

            if graph is not None:
                for parent_id, child_id in edges:
                    graph.addEdge(alias[parent_id], alias[child_id])
        finally:
            indicator.finish()

        self.graphId = 0
        self._graph_keys = list(self.graphDict.keys())

    @staticmethod
    def _parse_attributes(field: str, line_no: int) -> dict[str, str]:
        attrs: dict[str, str] = {}
        for pair in field.split(";"):
            if not pair:
                continue
            key, sep, value = pair.partition("=")
            if not sep or not key:
                raise ValueError(
                    f"Illegal attribute field '{field}' at line {line_no} in GFF file."
                )
            attrs[key] = value
        return attrs


__all__ = [
    "AS_KEY",
    "CHILD_REC",
    "GENE_REC",
    "ID_ATTR",
    "KNOWN_ATTRS",
    "PARENT_ATTR",
    "PARENT_REC",
    "SOURCE_NAME",
    "SpliceGraph",
    "SpliceGraphNode",
    "SpliceGraphParser",
]
