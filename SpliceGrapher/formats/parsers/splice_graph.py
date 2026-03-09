"""Parser boundary for networkx-backed splice graph GFF input."""

from __future__ import annotations

from collections.abc import Iterator
from typing import TextIO

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.core.splice_graph import (
    AS_KEY,
    END_CODON_KEY,
    ID_ATTR,
    ISO_KEY,
    KNOWN_ATTRS,
    PARENT_ATTR,
    START_CODON_KEY,
    VALID_GENES,
    VALID_RECTYPES,
    SpliceGraph,
)
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.progress import ProgressIndicator


class SpliceGraphParser:
    """Parse splice graph GFF files into networkx-backed ``SpliceGraph`` objects."""

    def __init__(self, fileRef: str | TextIO, *, verbose: bool = False) -> None:
        self.verbose = verbose
        self.instream = ez_open(fileRef) if isinstance(fileRef, str) else fileRef
        if self.instream is None:
            raise ValueError("No input file stream given.")
        self.graphDict: dict[str, SpliceGraph] = {}
        self.load_from_file()

    def __iter__(self) -> Iterator[SpliceGraph]:
        return iter(self.graphDict.values())

    def __len__(self) -> int:
        return len(self.graphDict)

    def load_from_file(self) -> None:
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


__all__ = ["SpliceGraphParser"]
