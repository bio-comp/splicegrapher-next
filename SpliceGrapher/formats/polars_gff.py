"""Optional Polars-backed helpers for GFF parsing and analytics."""

from __future__ import annotations

import importlib
from collections.abc import Iterable, Iterator
from pathlib import Path
from typing import TYPE_CHECKING

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType, Strand

if TYPE_CHECKING:
    import polars as pl

    from SpliceGrapher.formats.gene_model import GeneModel

_GFF_COLUMN_COUNT = 9
_GFF_COLUMNS = (
    "chrom",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes",
)

_RECORD_TYPE_BY_CASEFOLD = {record_type.value.casefold(): record_type for record_type in RecordType}


class PolarsNotInstalledError(ImportError):
    """Raised when optional Polars dependency is unavailable."""


def parse_gff_attributes(attribute_string: str) -> dict[str, str]:
    """Parse a GFF3 attribute column into key-value pairs.

    Malformed tokens (missing "=" or empty key/value) are ignored.
    """
    cleaned = attribute_string.strip()
    if not cleaned or cleaned == ".":
        return {}

    result: dict[str, str] = {}
    for token in cleaned.split(";"):
        item = token.strip()
        if not item or "=" not in item:
            continue

        key, value = item.split("=", 1)
        if not key or not value:
            continue

        result[key] = value

    return result


def _row_from_parts(parts: list[str]) -> dict[str, str | int | None]:
    if len(parts) < _GFF_COLUMN_COUNT:
        raise ValueError(f"expected {_GFF_COLUMN_COUNT} columns, got {len(parts)}")

    try:
        start = int(parts[3])
        end = int(parts[4])
    except ValueError as exc:
        raise ValueError("start/end must be integers") from exc

    attributes = parse_gff_attributes(parts[8])
    record_type = _coerce_record_type(parts[2])
    strand = coerce_enum(parts[6], Strand, field="strand")

    return {
        "chrom": parts[0],
        "source": parts[1],
        "type": record_type,
        "start": start,
        "end": end,
        "score": parts[5],
        "strand": strand,
        "phase": parts[7],
        "attributes": parts[8],
        "feature_id": attributes.get("ID"),
        "parent_id": attributes.get("Parent"),
        "name": attributes.get("Name"),
    }


def _coerce_record_type(value: str) -> RecordType:
    normalized = value.strip()
    try:
        return _RECORD_TYPE_BY_CASEFOLD[normalized.casefold()]
    except KeyError:
        return coerce_enum(normalized, RecordType, field="record_type")


def iter_gff_records(
    lines: Iterable[str], *, ignore_malformed: bool = False
) -> Iterator[dict[str, str | int | None]]:
    """Yield normalized row dictionaries from raw GFF lines."""
    for line_no, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue

        parts = stripped.split("\t")
        try:
            yield _row_from_parts(parts)
        except ValueError as exc:
            if ignore_malformed:
                continue
            raise ValueError(f"line {line_no}: {exc}") from exc


def load_gff_rows(
    path: str | Path,
    *,
    ignore_malformed: bool = False,
) -> list[dict[str, str | int | None]]:
    """Load GFF rows into dictionaries using the built-in parser path."""
    gff_path = Path(path)
    with gff_path.open("r", encoding="utf-8") as handle:
        return list(iter_gff_records(handle, ignore_malformed=ignore_malformed))


def load_gff_to_polars(
    path: str | Path,
    *,
    ignore_malformed: bool = False,
) -> "pl.DataFrame":
    """Load GFF rows into a Polars DataFrame (optional dependency)."""
    try:
        pl = importlib.import_module("polars")
    except ModuleNotFoundError as exc:
        raise PolarsNotInstalledError(
            "polars is not installed; install optional dependency to use this helper"
        ) from exc

    rows = load_gff_rows(path, ignore_malformed=ignore_malformed)
    schema = {
        "chrom": pl.Utf8,
        "source": pl.Utf8,
        "type": pl.Utf8,
        "start": pl.Int64,
        "end": pl.Int64,
        "score": pl.Utf8,
        "strand": pl.Utf8,
        "phase": pl.Utf8,
        "attributes": pl.Utf8,
        "feature_id": pl.Utf8,
        "parent_id": pl.Utf8,
        "name": pl.Utf8,
    }
    if not rows:
        return pl.DataFrame({name: [] for name in schema}, schema=schema)

    return pl.DataFrame(rows, schema=schema)


def _record_type_name(value: object, *, fallback: str) -> str:
    if isinstance(value, RecordType):
        return value.name
    value_string = str(value).strip()
    if not value_string:
        return fallback
    return value_string.upper()


def _feature_coordinates(feature: object) -> tuple[str, str, int, int]:
    """Extract ``(chromosome, strand, start, end)`` from legacy or composed models."""
    locus = getattr(feature, "locus", None)
    if locus is not None:
        return (
            str(getattr(locus, "chromosome")),
            str(getattr(locus, "strand")),
            int(getattr(locus, "minpos")),
            int(getattr(locus, "maxpos")),
        )
    return (
        str(getattr(feature, "chromosome")),
        str(getattr(feature, "strand")),
        int(getattr(feature, "minpos")),
        int(getattr(feature, "maxpos")),
    )


def _iter_gene_transcripts(gene: object) -> Iterator[object]:
    transcripts = getattr(gene, "transcripts", None)
    if isinstance(transcripts, dict):
        yield from transcripts.values()
        return

    for attr_name in ("mrna", "isoforms"):
        forms = getattr(gene, attr_name, None)
        if not isinstance(forms, dict):
            continue
        yield from forms.values()


def _iter_flattened_features(model: "GeneModel") -> Iterator[dict[str, str | int | None]]:
    """Yield flat feature rows from the gene model hierarchy."""
    for gene in model.iter_all_genes():
        gene_chrom, gene_strand, gene_start, gene_end = _feature_coordinates(gene)
        gene_attrs = getattr(gene, "attributes", {})
        gene_biotype = gene_attrs.get("gene_biotype") if isinstance(gene_attrs, dict) else None
        yield {
            "gene_id": str(getattr(gene, "id")),
            "transcript_id": None,
            "feature_type": "GENE",
            "chromosome": gene_chrom,
            "strand": gene_strand,
            "start_pos": gene_start,
            "end_pos": gene_end,
            "biotype": gene_biotype,
        }

        for transcript in _iter_gene_transcripts(gene):
            transcript_chrom, transcript_strand, transcript_start, transcript_end = (
                _feature_coordinates(transcript)
            )
            transcript_attrs = getattr(transcript, "attributes", {})
            transcript_biotype = (
                transcript_attrs.get("transcript_biotype")
                if isinstance(transcript_attrs, dict)
                else None
            )
            transcript_id = str(getattr(transcript, "id"))
            transcript_type = _record_type_name(
                getattr(transcript, "record_type", getattr(transcript, "feature_type", "")),
                fallback="TRANSCRIPT",
            )
            yield {
                "gene_id": str(getattr(gene, "id")),
                "transcript_id": transcript_id,
                "feature_type": transcript_type,
                "chromosome": transcript_chrom,
                "strand": transcript_strand,
                "start_pos": transcript_start,
                "end_pos": transcript_end,
                "biotype": transcript_biotype,
            }

            for exon in getattr(transcript, "exons", []):
                exon_chrom, exon_strand, exon_start, exon_end = _feature_coordinates(exon)
                yield {
                    "gene_id": str(getattr(gene, "id")),
                    "transcript_id": transcript_id,
                    "feature_type": "EXON",
                    "chromosome": exon_chrom,
                    "strand": exon_strand,
                    "start_pos": exon_start,
                    "end_pos": exon_end,
                    "biotype": None,
                }

            cds_regions = getattr(transcript, "cds_regions", getattr(transcript, "cds", []))
            for cds in cds_regions:
                cds_chrom, cds_strand, cds_start, cds_end = _feature_coordinates(cds)
                cds_type = _record_type_name(
                    getattr(cds, "record_type", getattr(cds, "feature_type", "")),
                    fallback="CDS",
                )
                yield {
                    "gene_id": str(getattr(gene, "id")),
                    "transcript_id": transcript_id,
                    "feature_type": cds_type,
                    "chromosome": cds_chrom,
                    "strand": cds_strand,
                    "start_pos": cds_start,
                    "end_pos": cds_end,
                    "biotype": None,
                }


def extract_to_dataframe(model: "GeneModel") -> "pl.DataFrame":
    """Build a typed Polars DataFrame from a flattened gene-model stream."""
    try:
        pl = importlib.import_module("polars")
    except ModuleNotFoundError as exc:
        raise PolarsNotInstalledError(
            "polars is not installed; install optional dependency to use this helper"
        ) from exc

    schema = {
        "gene_id": pl.Utf8,
        "transcript_id": pl.Utf8,
        "feature_type": pl.Utf8,
        "chromosome": pl.Categorical,
        "strand": pl.Categorical,
        "start_pos": pl.Int64,
        "end_pos": pl.Int64,
        "biotype": pl.Utf8,
    }
    rows = _iter_flattened_features(model)
    return pl.DataFrame(rows, schema=schema)


__all__ = [
    "PolarsNotInstalledError",
    "extract_to_dataframe",
    "iter_gff_records",
    "load_gff_rows",
    "load_gff_to_polars",
    "parse_gff_attributes",
]
