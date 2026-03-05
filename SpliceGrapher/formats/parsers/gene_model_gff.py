"""GFF parser boundary for ``GeneModel`` loading."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping, Sequence
from dataclasses import dataclass, field
from typing import TypeAlias

import structlog

# Safe import boundary: GeneModel loads this parser lazily at call time.
import SpliceGrapher.formats.gene_model as gm
from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.format_utils import comma_format
from SpliceGrapher.shared.progress import ProgressIndicator

RecordHandler: TypeAlias = Callable[["ParseContext", "ParsedRecord"], None]
LOGGER = structlog.get_logger(__name__)


@dataclass(slots=True)
class ParseStats:
    gene_count: int = 0
    exon_count: int = 0
    iso_count: int = 0
    mrna_count: int = 0
    cds_count: int = 0
    orphan_mrna_count: int = 0
    bad_lines: int = 0


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


@dataclass(slots=True)
class ParseContext:
    model: gm.GeneModel
    require_notes: bool
    verbose: bool
    ignore_errors: bool
    chromosomes: set[str] | None = None
    stats: ParseStats = field(default_factory=ParseStats)
    gene_alias: dict[str, str] = field(default_factory=dict)
    parent_cache: dict[tuple[str, str, bool, bool], gm.Gene | gm.Mrna] = field(default_factory=dict)
    parent_candidates: dict[str, tuple[str, ...]] = field(default_factory=dict)

    def fail(self, message: str) -> None:
        if self.ignore_errors:
            return
        raise RuntimeError(message)

    def write_verbose(self, message: str) -> None:
        if self.verbose:
            LOGGER.info("gene_model_gff_verbose", message=message)

    def report_bad_line(self, line_no: int) -> None:
        self.stats.bad_lines += 1
        self.write_verbose(
            f"line {line_no}: invalid GFF format (not enough columns); file may be corrupt"
        )
        if self.stats.bad_lines >= gm.MAX_BAD_LINES:
            LOGGER.error("gene_model_gff_invalid_input", line=line_no)
            raise ValueError("Invalid GFF input file")


def _subnames(full_string: str, delimiter: str) -> list[str]:
    parts = full_string.split(delimiter)
    return [delimiter.join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def _annotation_value(annots: Mapping[str, str], key: object) -> str | None:
    return annots.get(str(key))


def _normalize_chromosomes(chromosomes: Sequence[str] | str | None) -> set[str] | None:
    if chromosomes is None:
        return None
    if isinstance(chromosomes, str):
        return {chromosomes.lower()} if chromosomes else None
    return {str(chrom).lower() for chrom in chromosomes}


def _candidate_parent_ids(ctx: ParseContext, parent_string: str) -> tuple[str, ...]:
    cached = ctx.parent_candidates.get(parent_string)
    if cached is not None:
        return cached

    candidates: list[str] = []
    seen: set[str] = set()

    def add_candidate(candidate: str) -> None:
        if candidate and candidate not in seen:
            seen.add(candidate)
            candidates.append(candidate)

    add_candidate(parent_string)
    delimiters = [delimiter for delimiter in gm.FORM_DELIMITERS if delimiter in parent_string]
    for delimiter in delimiters:
        for candidate in _subnames(parent_string, delimiter):
            add_candidate(candidate)
    for item in parent_string.split(","):
        add_candidate(item)
        for piece in item.split("."):
            add_candidate(piece)

    result = tuple(candidates)
    ctx.parent_candidates[parent_string] = result
    return result


def _iter_record_lines(gff_records: gm.GffRecordSource, *, verbose: bool) -> Iterable[str]:
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


def _resolve_parent(
    ctx: ParseContext,
    parent_id: str,
    chrom: str,
    *,
    search_genes: bool = True,
    search_mrna: bool = True,
) -> gm.Gene | gm.Mrna | None:
    parent_string = parent_id.upper()
    chrom_key = chrom.lower()
    cache_key = (chrom_key, parent_string, search_genes, search_mrna)
    cached = ctx.parent_cache.get(cache_key)
    if cached is not None:
        return cached

    exact = ctx.model.get_parent(
        parent_string,
        chrom_key,
        search_genes=search_genes,
        search_mrna=search_mrna,
    )
    if exact is not None:
        ctx.parent_cache[cache_key] = exact
        return exact

    candidates = _candidate_parent_ids(ctx, parent_string)

    if search_mrna and chrom_key in ctx.model.mrna_forms:
        for candidate in candidates:
            transcript = ctx.model.mrna_forms[chrom_key].get(candidate)
            if transcript is not None:
                ctx.parent_cache[cache_key] = transcript
                return transcript

    if search_genes and chrom_key in ctx.model.model:
        for candidate in candidates:
            gene = ctx.model.model[chrom_key].get(candidate)
            if gene is not None:
                ctx.parent_cache[cache_key] = gene
                return gene
    return None


def _resolve_strand(
    ctx: ParseContext,
    *,
    observed: str,
    expected: str,
    mismatch_message: str,
) -> str:
    if observed in gm.VALID_STRANDS and observed != expected:
        ctx.fail(mismatch_message)
        return observed
    return expected


def _parse_record_line(ctx: ParseContext, line: str, line_no: int) -> ParsedRecord | None:
    stripped = line.rstrip()
    if not stripped or stripped[0] == "#":
        return None

    parts = stripped.split("\t")
    if len(parts) < 7:
        ctx.report_bad_line(line_no)
        return None

    annots = ctx.model.get_annotation_dict(parts[-1])
    chrom = parts[0].lower()
    if ctx.chromosomes and chrom not in ctx.chromosomes:
        return None

    try:
        rec_type = coerce_enum(parts[2].lower(), RecordType, field="record_type")
    except ValueError as exc:
        raise ValueError(f"line {line_no}: unknown record type '{parts[2]}'") from exc
    rec_type = gm.RECTYPE_MAP.get(rec_type, rec_type)

    start_pos = int(parts[3])
    end_pos = int(parts[4])
    try:
        strand = coerce_enum(parts[6], Strand, field="strand").value
    except ValueError as exc:
        raise ValueError(f"line {line_no}: unknown strand '{parts[6]}'") from exc

    ctx.model.found_types[rec_type] = (
        rec_type in gm.KNOWN_RECTYPES and rec_type not in gm.IGNORE_RECTYPES
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


def _handle_gene_record(ctx: ParseContext, record: ParsedRecord) -> None:
    gid = _annotation_value(record.annots, gm.ID_FIELD)
    if gid is None:
        raise ValueError(
            f"line {record.line_no}: {record.rec_type} record has no ID field:\n{record.raw_line}\n"
        )
    gid = gid.upper()

    name = ctx.model.get_annotation(gm.NAME_FIELD, record.annots, None)
    if name:
        name = ctx.model.clean_name(name)

    note = ctx.model.get_annotation(gm.NOTE_FIELD, record.annots)
    if not note and ctx.require_notes:
        return

    if record.strand not in gm.VALID_STRANDS:
        ctx.fail(f"line {record.line_no}: {record.rec_type} record with unknown strand")

    if record.chrom not in ctx.model.model:
        ctx.model.model[record.chrom] = {}
        ctx.model.add_chromosome(1, record.end_pos, record.chrom)

    if record.rec_type == RecordType.PSEUDOGENE:
        gene_obj: gm.Gene = gm.PseudoGene(
            gid,
            note,
            record.start_pos,
            record.end_pos,
            record.chrom,
            record.strand,
            name,
            record.annots,
        )
    else:
        gene_obj = gm.Gene(
            gid,
            note,
            record.start_pos,
            record.end_pos,
            record.chrom,
            record.strand,
            name,
            record.annots,
        )

    other = ctx.model.all_genes.get(str(gene_obj))
    if other is not None:
        ctx.fail(
            f"line {record.line_no}: gene {gene_obj.id} associated with multiple loci: "
            f"{other.minpos}-{other.maxpos} and {record.start_pos}-{record.end_pos}"
        )

    ctx.model.all_chr[record.chrom].update(gene_obj)
    ctx.model.add_gene(gene_obj)
    ctx.gene_alias[gene_obj.name.upper()] = gene_obj.id.upper()
    ctx.stats.gene_count += 1


def _resolve_exon_parent(ctx: ParseContext, record: ParsedRecord) -> gm.Gene | None:
    parent_record: gm.Gene | gm.Mrna | None = None
    tried: set[str] = set()
    for key in gm.POSSIBLE_GENE_FIELDS:
        parent_name = _annotation_value(record.annots, key)
        if parent_name is None or parent_name in tried:
            continue
        parent_record = _resolve_parent(ctx, parent_name, record.chrom)
        if parent_record is not None:
            break
        tried.add(parent_name)

    if parent_record is None:
        return None
    if isinstance(parent_record, gm.Mrna):
        parent_gene = ctx.model.get_mrna_parent(record.chrom, parent_record.id)
        if parent_gene is None and isinstance(parent_record.parent, gm.Gene):
            parent_gene = parent_record.parent
        if parent_gene is None:
            ctx.fail(f"line {record.line_no}: Mrna parent is missing gene for exon record")
            return None
        return parent_gene
    if not isinstance(parent_record, gm.Gene):
        ctx.fail(f"line {record.line_no}: exon parent {parent_record} is not a gene")
        return None
    return parent_record


def _resolve_isoform(
    ctx: ParseContext,
    record: ParsedRecord,
    gene_obj: gm.Gene,
) -> tuple[gm.Isoform | None, bool]:
    isoform = None
    iso_name = ""
    tried: set[str] = set()
    for key in gm.POSSIBLE_FORM_FIELDS:
        name = _annotation_value(record.annots, key)
        if name is None or name in tried:
            continue
        iso_name = name
        if iso_name in gene_obj.isoforms:
            isoform = gene_obj.isoforms[iso_name]
            break
        tried.add(iso_name)

    if not (isoform or iso_name):
        return (None, False)
    if isoform is not None:
        return (isoform, False)

    strand = _resolve_strand(
        ctx,
        observed=record.strand,
        expected=gene_obj.strand,
        mismatch_message=(
            f"line {record.line_no}: exon strand ({record.strand}) != gene "
            f"strand ({gene_obj.strand}) for {gene_obj.id}"
        ),
    )
    iso_attr = {
        str(gm.PARENT_FIELD): gene_obj.id,
        str(gm.NAME_FIELD): iso_name,
        str(gm.ID_FIELD): iso_name,
    }
    return (
        gm.Isoform(
            iso_name,
            record.start_pos,
            record.end_pos,
            record.chrom,
            strand,
            attr=iso_attr,
        ),
        True,
    )


def _handle_exon_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if record.chrom not in ctx.model.model:
        return

    gene_obj = _resolve_exon_parent(ctx, record)
    if gene_obj is None:
        return

    isoform, created_isoform = _resolve_isoform(ctx, record, gene_obj)
    if isoform is None:
        return

    strand = _resolve_strand(
        ctx,
        observed=record.strand,
        expected=gene_obj.strand,
        mismatch_message=(
            f"line {record.line_no}: exon strand ({record.strand}) != gene "
            f"strand ({gene_obj.strand}) for {gene_obj.id}"
        ),
    )
    exon = gm.Exon(record.start_pos, record.end_pos, record.chrom, strand, record.annots)
    if gene_obj.add_exon(isoform, exon):
        ctx.stats.exon_count += 1
    if created_isoform:
        ctx.stats.iso_count += 1


def _resolve_mrna_gene(ctx: ParseContext, record: ParsedRecord, parent_id: str) -> gm.Gene | None:
    parent_id_upper = parent_id.upper()
    parent_candidate = _resolve_parent(ctx, parent_id_upper, record.chrom)
    if isinstance(parent_candidate, gm.Gene):
        return parent_candidate

    alias = ctx.gene_alias.get(parent_id_upper)
    if alias is None:
        return None
    alias_candidate = _resolve_parent(ctx, alias, record.chrom)
    if isinstance(alias_candidate, gm.Gene):
        return alias_candidate
    return None


def _handle_mrna_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if record.chrom not in ctx.model.model:
        ctx.fail(
            f"line {record.line_no}: Mrna with missing chromosome dictionary {record.chrom} "
            f"(known: {','.join(ctx.model.model.keys())})"
        )
        return

    transcript_id = _annotation_value(record.annots, gm.ID_FIELD)
    if transcript_id is None:
        ctx.fail(f"line {record.line_no}: Mrna with missing ID")
        return
    transcript_id = transcript_id.upper()

    parent_id = ctx.model.get_annotation(gm.PARENT_FIELD, record.annots)
    if parent_id is None:
        return

    mrna_gene = _resolve_mrna_gene(ctx, record, parent_id)
    if mrna_gene is None:
        ctx.stats.orphan_mrna_count += 1
        parent_id_upper = parent_id.upper()
        alias = ctx.gene_alias.get(parent_id_upper, parent_id_upper)
        if ctx.verbose:
            if alias == parent_id_upper:
                ctx.write_verbose(
                    f"line {record.line_no}: no gene {parent_id_upper} found for {record.rec_type}"
                )
            else:
                ctx.write_verbose(
                    f"line {record.line_no}: no gene '{parent_id_upper}' or '{alias}' "
                    f"found for {record.rec_type}"
                )
        ctx.fail(f"line {record.line_no}: no gene parent found for transcript {parent_id_upper}")
        return

    strand = _resolve_strand(
        ctx,
        observed=record.strand,
        expected=mrna_gene.strand,
        mismatch_message=(
            f"line {record.line_no}: Mrna strand ({record.strand}) does not "
            f"match gene strand ({mrna_gene.strand})"
        ),
    )
    mrna_attr = {
        str(gm.PARENT_FIELD): mrna_gene.id,
        str(gm.NAME_FIELD): transcript_id,
        str(gm.ID_FIELD): transcript_id,
    }
    mrna = gm.Mrna(
        transcript_id,
        record.start_pos,
        record.end_pos,
        record.chrom,
        strand,
        attr=mrna_attr,
    )
    mrna_gene.add_mrna(mrna)
    ctx.model.mrna_forms.setdefault(record.chrom, {})[transcript_id] = mrna
    ctx.model.mrna_gene.setdefault(record.chrom, {})[transcript_id] = mrna_gene
    ctx.stats.mrna_count += 1


def _handle_transcript_region_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if record.chrom not in ctx.model.model:
        ctx.fail(
            f"line {record.line_no}: {record.rec_type} has unrecognized chromosome: "
            f"{record.chrom} (known: {','.join(ctx.model.model.keys())})"
        )
        return
    if record.chrom not in ctx.model.mrna_forms:
        ctx.fail(
            f"line {record.line_no}: {record.rec_type} has unrecognized chromosome: "
            f"{record.chrom} (known: {','.join(ctx.model.mrna_forms.keys())})"
        )
        return

    parent_id = _annotation_value(record.annots, gm.PARENT_FIELD)
    if parent_id is None:
        return
    mrna_record = _resolve_parent(ctx, parent_id, record.chrom, search_genes=False)
    if mrna_record is None:
        ctx.write_verbose(f"line {record.line_no}: no Mrna {parent_id} found for {record.rec_type}")
        return
    if not isinstance(mrna_record, gm.Mrna):
        ctx.fail(f"line {record.line_no}: parent {parent_id} is not an Mrna record")
        return
    parent_gene = ctx.model.get_mrna_parent(record.chrom, mrna_record.id)
    if parent_gene is None and isinstance(mrna_record.parent, gm.Gene):
        parent_gene = mrna_record.parent
    if parent_gene is None:
        ctx.fail(
            f"line {record.line_no}: Mrna {mrna_record.id} is missing a parent gene for CDS record"
        )
        return

    strand = _resolve_strand(
        ctx,
        observed=record.strand,
        expected=mrna_record.strand,
        mismatch_message=(
            f"line {record.line_no}: CDS strand ({record.strand}) does not "
            f"match Mrna strand ({mrna_record.strand})"
        ),
    )
    cds = gm.cds_factory(
        record.rec_type,
        record.start_pos,
        record.end_pos,
        record.chrom,
        strand,
        record.annots,
    )
    if parent_gene.add_cds(mrna_record, cds):
        ctx.stats.cds_count += 1


def _extract_parent_gene_id(parent_id: str) -> str | None:
    if not parent_id:
        return None
    primary_parent = parent_id.split(",", 1)[0].strip()
    if not primary_parent:
        return None
    return primary_parent.split(".", 1)[0]


def _handle_misc_feature_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if record.chrom not in ctx.model.model:
        return
    parent_value = _annotation_value(record.annots, gm.PARENT_FIELD)
    if parent_value is None:
        return
    parent_gene_id = _extract_parent_gene_id(parent_value)
    if not parent_gene_id:
        return
    feature_gene = ctx.model.model[record.chrom].get(parent_gene_id)
    if feature_gene is None:
        return
    try:
        feature_gene.add_feature(
            gm.BaseFeature(
                record.rec_type,
                record.start_pos,
                record.end_pos,
                record.chrom,
                record.strand,
                record.annots,
            )
        )
    except ValueError as err:
        ctx.fail(f"line {record.line_no}: {err}")


def _handle_chromosome_record(ctx: ParseContext, record: ParsedRecord) -> None:
    ctx.model.add_chromosome(record.start_pos, record.end_pos, record.chrom)


RECORD_HANDLERS: dict[RecordType, RecordHandler] = {
    RecordType.GENE: _handle_gene_record,
    RecordType.PSEUDOGENE: _handle_gene_record,
    RecordType.EXON: _handle_exon_record,
    RecordType.PSEUDOGENIC_EXON: _handle_exon_record,
    RecordType.MRNA: _handle_mrna_record,
    RecordType.PSEUDOGENIC_TRANSCRIPT: _handle_mrna_record,
    RecordType.CDS: _handle_transcript_region_record,
    RecordType.FIVE_PRIME_UTR: _handle_transcript_region_record,
    RecordType.THREE_PRIME_UTR: _handle_transcript_region_record,
    RecordType.CHROMOSOME: _handle_chromosome_record,
}


def load_gene_model_records(
    model: gm.GeneModel,
    gff_records: gm.GffRecordSource,
    *,
    require_notes: bool = False,
    chromosomes: Sequence[str] | str | None = None,
    verbose: bool = False,
    ignore_errors: bool = False,
) -> None:
    """Load gene model records from GFF-like input into ``model``."""
    chrom_filter = _normalize_chromosomes(chromosomes)
    if verbose and chrom_filter is not None:
        LOGGER.info(
            "gene_model_gff_chromosome_filter",
            chromosomes=sorted(chrom_filter),
        )

    model.model = {}
    model.mrna_forms = {}
    model.mrna_gene = {}
    model.all_genes = {}
    model.found_types = {}
    model.all_chr = {}
    model.chromosome_index = {}
    ctx = ParseContext(
        model=model,
        require_notes=require_notes,
        verbose=verbose,
        ignore_errors=ignore_errors,
        chromosomes=chrom_filter,
    )

    indicator = ProgressIndicator(1000000, verbose=verbose)
    for line_no, line in enumerate(_iter_record_lines(gff_records, verbose=verbose), start=1):
        indicator.update()
        record = _parse_record_line(ctx, line, line_no)
        if record is None:
            continue
        handler = RECORD_HANDLERS.get(record.rec_type, _handle_misc_feature_record)
        handler(ctx, record)

    indicator.finish()
    if verbose:
        if ctx.stats.gene_count > 0:
            LOGGER.info(
                "gene_model_gff_load_summary",
                gene_count=ctx.stats.gene_count,
                isoform_count=ctx.stats.iso_count,
                exon_count=ctx.stats.exon_count,
                exon_per_gene=round(float(ctx.stats.exon_count) / ctx.stats.gene_count, 1),
                mrna_count=ctx.stats.mrna_count,
                orphan_mrna_count=ctx.stats.orphan_mrna_count,
                cds_count=ctx.stats.cds_count,
                cds_per_gene=round(float(ctx.stats.cds_count / ctx.stats.gene_count), 1),
                gene_count_fmt=comma_format(ctx.stats.gene_count),
                isoform_count_fmt=comma_format(ctx.stats.iso_count),
                exon_count_fmt=comma_format(ctx.stats.exon_count),
                mrna_count_fmt=comma_format(ctx.stats.mrna_count),
                orphan_mrna_count_fmt=comma_format(ctx.stats.orphan_mrna_count),
                cds_count_fmt=comma_format(ctx.stats.cds_count),
            )
        else:
            LOGGER.warning("gene_model_gff_no_genes_loaded")


__all__ = ["load_gene_model_records"]
