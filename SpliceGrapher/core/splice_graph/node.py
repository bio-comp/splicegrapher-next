"""Node models for the splice-graph core."""

from __future__ import annotations

from dataclasses import dataclass, field

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AlternativeSplicingEvent, AlternativeSplicingEventName, Strand

from .constants import (
    AS_KEY,
    DISPOSITION_KEY,
    END_CODON_KEY,
    ISO_KEY,
    KNOWN_NODE,
    PREDICTED_NODE,
    PUTATIVE_CHILDREN,
    PUTATIVE_PARENTS,
    START_CODON_KEY,
    UNRESOLVED_NODE,
    NodeAttributeValue,
    coerce_alt_splicing_event,
)


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
    alt_form_set: set[AlternativeSplicingEvent]
    isoform_set: set[str]
    orig_start: int
    orig_end: int

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
        self.alt_form_set = set()
        self.isoform_set = set()
        self.orig_start = self.start
        self.orig_end = self.end

    def acceptor_end(self) -> int:
        return self.start

    def donor_end(self) -> int:
        return self.end

    def _sync_alt_form_attr(self) -> None:
        self.attrs[AS_KEY] = ",".join(sorted(form.value for form in self.alt_form_set))

    def add_alt_form(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        if isinstance(form, str) and not form.strip():
            return
        event = coerce_alt_splicing_event(form)
        self.alt_form_set.add(event)
        self._sync_alt_form_attr()

    def remove_alt_form(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        if isinstance(form, str) and not form.strip():
            return
        event = coerce_alt_splicing_event(form)
        self.alt_form_set.discard(event)
        self._sync_alt_form_attr()

    def add_attribute(self, key: str, value: NodeAttributeValue) -> None:
        self.attrs[key] = value

    def add_codon(self, codon: tuple[int, int], codon_type: str) -> None:
        if len(codon) != 2:
            raise ValueError(f"Codons must be 2-tuples; received {codon!r}")
        pos = min(codon) if self.strand == Strand.PLUS.value else max(codon)
        if self.contains(pos):
            existing = self.attrs.setdefault(codon_type, set())
            if not isinstance(existing, set):
                raise TypeError(f"Codon attribute {codon_type} is not a set")
            existing.add(pos)

    def add_start_codon(self, codon: tuple[int, int]) -> None:
        self.add_codon(codon, START_CODON_KEY)

    def add_end_codon(self, codon: tuple[int, int]) -> None:
        self.add_codon(codon, END_CODON_KEY)

    def add_forms_from_string(self, forms: str) -> None:
        for form in forms.split(","):
            self.add_alt_form(form.strip())

    def add_isoform(self, isoform: str) -> None:
        if isoform is None:
            raise ValueError(f"Received illegal isoform for {self.id}")
        self.isoform_set.add(isoform)
        self.attrs[ISO_KEY] = ",".join(sorted(self.isoform_set))

    def add_isoform_string(self, isoform_string: str) -> None:
        for isoform in isoform_string.split(","):
            cleaned = isoform.strip()
            if cleaned:
                self.add_isoform(cleaned)

    def alt_forms(self) -> list[str]:
        return [form.value for form in sorted(self.alt_form_set, key=lambda form: form.value)]

    def alt_form_string(self) -> str:
        value = self.attrs.get(AS_KEY)
        return value if isinstance(value, str) else ""

    def attribute_string(self) -> str:
        attr_parts: list[str] = []
        for key in sorted(self.attrs):
            if key == AS_KEY and not self.alt_form_set:
                continue
            if key == ISO_KEY and not self.isoform_set:
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

    def codon_string(self, codon_type: str) -> str | None:
        codons = self.codons(codon_type)
        if not codons:
            return None
        return ",".join(str(value) for value in codons)

    def contains(self, pos: int) -> bool:
        return self.minpos <= pos <= self.maxpos

    def downstream_of(self, pos: int) -> bool:
        return self.minpos > pos if self.strand == Strand.PLUS.value else self.maxpos < pos

    def end_codons(self) -> list[int]:
        return self.codons(END_CODON_KEY)

    def end_codon_string(self) -> str | None:
        return self.codon_string(END_CODON_KEY)

    def has_as(self) -> bool:
        return bool(self.alt_form_set)

    def has_disposition(self, disposition: str) -> bool:
        value = self.attrs.get(DISPOSITION_KEY)
        return value == disposition

    def is_alt_acceptor(self) -> bool:
        return AlternativeSplicingEvent.ALT3 in self.alt_form_set

    def is_alt_donor(self) -> bool:
        return AlternativeSplicingEvent.ALT5 in self.alt_form_set

    def is_known(self) -> bool:
        return self.has_disposition(KNOWN_NODE)

    def is_predicted(self) -> bool:
        return self.has_disposition(PREDICTED_NODE)

    def is_retained_intron(self) -> bool:
        return AlternativeSplicingEvent.IR in self.alt_form_set

    def is_skipped_exon(self) -> bool:
        return AlternativeSplicingEvent.ES in self.alt_form_set

    def is_unresolved(self) -> bool:
        return self.has_disposition(UNRESOLVED_NODE)

    def isoform_list(self) -> list[str]:
        return sorted(self.isoform_set)

    def isoform_string(self) -> str | None:
        value = self.attrs.get(ISO_KEY)
        return value if isinstance(value, str) else None

    def putative_children(self) -> set[str]:
        value = self.attrs.get(PUTATIVE_CHILDREN)
        if isinstance(value, str) and value:
            return set(value.split(","))
        return set()

    def putative_parents(self) -> set[str]:
        value = self.attrs.get(PUTATIVE_PARENTS)
        if isinstance(value, str) and value:
            return set(value.split(","))
        return set()

    def start_codons(self) -> list[int]:
        return self.codons(START_CODON_KEY)

    def start_codon_string(self) -> str | None:
        return self.codon_string(START_CODON_KEY)

    def update(self, minpos: int, maxpos: int) -> None:
        self.minpos = minpos
        self.maxpos = maxpos
        if self.strand == Strand.MINUS.value:
            self.start = maxpos
            self.end = minpos
        else:
            self.start = minpos
            self.end = maxpos

    def upstream_of(self, pos: int) -> bool:
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


__all__ = ["NullNode", "SpliceGraphNode"]
