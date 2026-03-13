"""Record parsing helpers for ``GeneModel`` GFF loading."""

from __future__ import annotations

from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from urllib.parse import unquote

import structlog

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.shared.file_utils import ez_open

LOGGER = structlog.get_logger(__name__)


@dataclass(slots=True)
class ParsedRecord:
    line_no: int
    raw_line: str
    annots: dict[str, str]
    chrom: str
    rec_type: RecordType
    start_pos: int
    end_pos: int
    strand: str


def annotation_value(annots: Mapping[str, str], key: str) -> str | None:
    return annots.get(key)


def normalize_chromosomes(chromosomes: Sequence[str] | str | None) -> set[str] | None:
    if chromosomes is None:
        return None
    if isinstance(chromosomes, str):
        return {chromosomes.lower()} if chromosomes else None
    return {str(chrom).lower() for chrom in chromosomes}


def clean_name(value: str) -> str:
    revised = unquote(value)
    revised = revised.replace(",", "")
    return revised.replace(" ", "-")


def parent_annotation(record: ParsedRecord) -> str | None:
    return annotation_value(record.annots, model_domain.PARENT_FIELD)


def parse_annotations(annotation_string: str) -> dict[str, str]:
    value_string = annotation_string.replace(" ", "")
    result: dict[str, str] = {}
    for item in value_string.split(";"):
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        if not key:
            continue
        result[key] = value
    return result


def iter_record_lines(
    gff_records: model_domain.GffRecordSource,
    *,
    verbose: bool,
) -> Iterable[str]:
    if isinstance(gff_records, str):
        if verbose:
            LOGGER.info("gene_model_gff_loading_file", source=gff_records)
        return ez_open(gff_records)
    if isinstance(gff_records, (list, set, tuple)):
        if verbose:
            LOGGER.info(
                "gene_model_gff_loading_records",
                record_count=len(gff_records),
            )
        return gff_records
    raise ValueError(
        "Unrecognized GFF record source "
        f"({type(gff_records).__name__}); must be file path or a list/set of strings."
    )


def parse_record_line(ctx: ParseContext, line: str, line_no: int) -> ParsedRecord | None:
    stripped = line.rstrip()
    if not stripped or stripped[0] == "#":
        return None

    parts = stripped.split("\t")
    if len(parts) < 7:
        ctx.report_bad_line(line_no)
        return None

    annots = parse_annotations(parts[-1])
    chrom = parts[0].lower()
    if ctx.chromosomes and chrom not in ctx.chromosomes:
        return None

    try:
        rec_type = coerce_enum(parts[2].lower(), RecordType, field="record_type")
    except ValueError as exc:
        raise ValueError(f"line {line_no}: unknown record type '{parts[2]}'") from exc
    rec_type = model_domain.RECTYPE_MAP.get(rec_type, rec_type)

    start_pos = int(parts[3])
    end_pos = int(parts[4])
    try:
        strand = coerce_enum(parts[6], Strand, field="strand").value
    except ValueError as exc:
        raise ValueError(f"line {line_no}: unknown strand '{parts[6]}'") from exc

    ctx.model.found_types[rec_type] = (
        rec_type in model_domain.KNOWN_RECTYPES and rec_type not in model_domain.IGNORE_RECTYPES
    )
    if not ctx.model.found_types[rec_type]:
        return None

    return ParsedRecord(
        line_no=line_no,
        raw_line=line,
        annots=annots,
        chrom=chrom,
        rec_type=rec_type,
        start_pos=start_pos,
        end_pos=end_pos,
        strand=strand,
    )


__all__ = [
    "ParsedRecord",
    "annotation_value",
    "clean_name",
    "iter_record_lines",
    "normalize_chromosomes",
    "parent_annotation",
    "parse_annotations",
    "parse_record_line",
]
