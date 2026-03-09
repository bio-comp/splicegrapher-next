"""Record handlers for ``GeneModel`` GFF loading."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from typing import TypeAlias

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_records import (
    ParsedRecord,
    annotation_value,
    clean_name,
    parent_annotation,
)

RecordHandler: TypeAlias = Callable[[ParseContext, ParsedRecord], None]


def _subnames(full_string: str, delimiter: str) -> list[str]:
    parts = full_string.split(delimiter)
    return [delimiter.join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def _known_chromosomes(mapping: Mapping[str, object]) -> str:
    return ",".join(mapping.keys())


def _has_chromosome(
    ctx: ParseContext,
    record: ParsedRecord,
    mapping: Mapping[str, object],
    *,
    fail_message: str | None = None,
) -> bool:
    if record.chrom in mapping:
        return True
    if fail_message is not None:
        ctx.fail(fail_message)
    return False


def _pending_key(chrom: str, transcript_id: str) -> tuple[str, str]:
    return (chrom, transcript_id.upper())


def _queue_pending_child(
    pending_map: dict[tuple[str, str], list[ParsedRecord]],
    record: ParsedRecord,
    parent_id: str,
) -> None:
    pending_map.setdefault(_pending_key(record.chrom, parent_id), []).append(record)


def _drain_pending_children(
    ctx: ParseContext,
    *,
    transcript: model_domain.Transcript,
    parent_gene: model_domain.Gene,
    chrom: str,
) -> None:
    key = _pending_key(chrom, transcript.id)

    for pending in ctx.pending_exons.pop(key, []):
        strand = resolve_strand(
            ctx,
            observed=pending.strand,
            expected=transcript.strand,
            mismatch_message=(
                f"line {pending.line_no}: exon strand ({pending.strand}) != transcript "
                f"strand ({transcript.strand}) for {transcript.id}"
            ),
        )
        exon = model_domain.Exon(
            pending.start_pos,
            pending.end_pos,
            pending.chrom,
            strand,
            pending.annots,
        )
        if parent_gene.add_exon(transcript, exon):
            ctx.stats.exon_count += 1

    for pending in ctx.pending_regions.pop(key, []):
        strand = resolve_strand(
            ctx,
            observed=pending.strand,
            expected=transcript.strand,
            mismatch_message=(
                f"line {pending.line_no}: {pending.rec_type} strand ({pending.strand}) != "
                f"transcript strand ({transcript.strand}) for {transcript.id}"
            ),
        )
        cds = model_domain.cds_factory(
            pending.rec_type,
            pending.start_pos,
            pending.end_pos,
            pending.chrom,
            strand,
            pending.annots,
        )
        if parent_gene.add_cds(transcript, cds):
            ctx.stats.cds_count += 1


def candidate_parent_ids(ctx: ParseContext, parent_string: str) -> tuple[str, ...]:
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
    delimiters = [
        delimiter for delimiter in model_domain.FORM_DELIMITERS if delimiter in parent_string
    ]
    for delimiter in delimiters:
        for candidate in _subnames(parent_string, delimiter):
            add_candidate(candidate)
    for item in parent_string.split(","):
        add_candidate(item)
        for piece in item.split("."):
            add_candidate(piece)

    result = tuple(candidates)
    ctx.cache_parent_candidate(parent_string, result)
    return result


def resolve_parent(
    ctx: ParseContext,
    parent_id: str,
    chrom: str,
    *,
    search_genes: bool = True,
    search_mrna: bool = True,
) -> model_domain.Gene | model_domain.Transcript | None:
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
        ctx.cache_parent(cache_key, exact)
        return exact

    candidates = candidate_parent_ids(ctx, parent_string)

    if search_mrna and chrom_key in ctx.model.mrna_forms:
        for candidate in candidates:
            transcript = ctx.model.mrna_forms[chrom_key].get(candidate)
            if transcript is not None:
                ctx.cache_parent(cache_key, transcript)
                return transcript

    if search_genes and chrom_key in ctx.model.model:
        for candidate in candidates:
            gene = ctx.model.model[chrom_key].get(candidate)
            if gene is not None:
                ctx.cache_parent(cache_key, gene)
                return gene
    return None


def resolve_strand(
    ctx: ParseContext,
    *,
    observed: str,
    expected: str,
    mismatch_message: str,
) -> str:
    if observed in model_domain.VALID_STRANDS and observed != expected:
        ctx.fail(mismatch_message)
        return observed
    return expected


def _resolve_exon_parent(ctx: ParseContext, record: ParsedRecord) -> model_domain.Gene | None:
    parent_record: model_domain.Gene | model_domain.Transcript | None = None
    tried: set[str] = set()
    for key in model_domain.POSSIBLE_GENE_FIELDS:
        parent_name = annotation_value(record.annots, key)
        if parent_name is None or parent_name in tried:
            continue
        parent_record = resolve_parent(ctx, parent_name, record.chrom)
        if parent_record is not None:
            break
        tried.add(parent_name)

    if parent_record is None:
        return None
    if isinstance(parent_record, model_domain.Transcript):
        parent_gene = ctx.model.get_mrna_parent(record.chrom, parent_record.id)
        if parent_gene is None:
            ctx.fail(f"line {record.line_no}: Mrna parent is missing gene for exon record")
            return None
        return parent_gene
    if not isinstance(parent_record, model_domain.Gene):
        ctx.fail(f"line {record.line_no}: exon parent {parent_record} is not a gene")
        return None
    return parent_record


def _resolve_isoform(
    ctx: ParseContext,
    record: ParsedRecord,
    gene_obj: model_domain.Gene,
) -> tuple[model_domain.Transcript | None, bool]:
    isoform = None
    iso_name = ""
    tried: set[str] = set()
    for key in model_domain.POSSIBLE_FORM_FIELDS:
        name = annotation_value(record.annots, key)
        if name is None or name in tried:
            continue
        iso_name = name
        if iso_name in gene_obj.transcripts:
            existing = gene_obj.transcripts[iso_name]
            if isinstance(existing, model_domain.Transcript):
                isoform = existing
            else:
                isoform = model_domain.Transcript(
                    iso_name,
                    existing.minpos,
                    existing.maxpos,
                    existing.chromosome,
                    existing.strand,
                    attr=existing.attributes,
                )
                for exon in existing.exons:
                    isoform.add_exon(exon)
                gene_obj.add_transcript(isoform)
            break
        tried.add(iso_name)

    if not (isoform or iso_name):
        return (None, False)
    if isoform is not None:
        return (isoform, False)

    strand = resolve_strand(
        ctx,
        observed=record.strand,
        expected=gene_obj.strand,
        mismatch_message=(
            f"line {record.line_no}: exon strand ({record.strand}) != gene "
            f"strand ({gene_obj.strand}) for {gene_obj.id}"
        ),
    )
    iso_attr = {
        str(model_domain.PARENT_FIELD): gene_obj.id,
        str(model_domain.NAME_FIELD): iso_name,
        str(model_domain.ID_FIELD): iso_name,
    }
    return (
        model_domain.Transcript(
            iso_name,
            record.start_pos,
            record.end_pos,
            record.chrom,
            strand,
            attr=iso_attr,
        ),
        True,
    )


def _resolve_mrna_gene(
    ctx: ParseContext,
    record: ParsedRecord,
    parent_id: str,
) -> model_domain.Gene | None:
    parent_id_upper = parent_id.upper()
    parent_candidate = resolve_parent(ctx, parent_id_upper, record.chrom)
    if isinstance(parent_candidate, model_domain.Gene):
        return parent_candidate

    alias = ctx.gene_alias.get(parent_id_upper)
    if alias is None:
        return None
    alias_candidate = resolve_parent(ctx, alias, record.chrom)
    if isinstance(alias_candidate, model_domain.Gene):
        return alias_candidate
    return None


def handle_gene_record(ctx: ParseContext, record: ParsedRecord) -> None:
    gid = annotation_value(record.annots, model_domain.ID_FIELD)
    if gid is None:
        raise ValueError(
            f"line {record.line_no}: {record.rec_type} record has no ID field:\n{record.raw_line}\n"
        )
    gid = gid.upper()

    name = annotation_value(record.annots, model_domain.NAME_FIELD)
    if name:
        name = clean_name(name)

    note = annotation_value(record.annots, model_domain.NOTE_FIELD)
    if not note and ctx.require_notes:
        return

    if record.strand not in model_domain.VALID_STRANDS:
        ctx.fail(f"line {record.line_no}: {record.rec_type} record with unknown strand")

    if record.chrom not in ctx.model.model:
        ctx.model.model[record.chrom] = {}
        ctx.model.add_chromosome(1, record.end_pos, record.chrom)

    if record.rec_type == RecordType.PSEUDOGENE:
        gene_obj: model_domain.Gene = model_domain.PseudoGene(
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
        gene_obj = model_domain.Gene(
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
    alias_key = (gene_obj.name or gene_obj.id).upper()
    ctx.gene_alias[alias_key] = gene_obj.id.upper()
    ctx.stats.gene_count += 1


def handle_exon_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not _has_chromosome(ctx, record, ctx.model.model):
        return

    parent_id = parent_annotation(record)
    if parent_id is None:
        return
    existing_transcript = resolve_parent(ctx, parent_id, record.chrom, search_genes=False)
    if isinstance(existing_transcript, model_domain.Transcript):
        parent_gene = ctx.model.get_mrna_parent(record.chrom, existing_transcript.id)
        if parent_gene is None:
            ctx.fail(
                f"line {record.line_no}: transcript {existing_transcript.id} is missing parent gene"
            )
            return
        strand = resolve_strand(
            ctx,
            observed=record.strand,
            expected=existing_transcript.strand,
            mismatch_message=(
                f"line {record.line_no}: exon strand ({record.strand}) != transcript "
                f"strand ({existing_transcript.strand}) for {existing_transcript.id}"
            ),
        )
        exon = model_domain.Exon(
            record.start_pos,
            record.end_pos,
            record.chrom,
            strand,
            record.annots,
        )
        if parent_gene.add_exon(existing_transcript, exon):
            ctx.stats.exon_count += 1
        return

    gene_obj = _resolve_exon_parent(ctx, record)
    if gene_obj is None:
        _queue_pending_child(ctx.pending_exons, record, parent_id)
        return

    isoform, created_isoform = _resolve_isoform(ctx, record, gene_obj)
    if isoform is None:
        return

    strand = resolve_strand(
        ctx,
        observed=record.strand,
        expected=gene_obj.strand,
        mismatch_message=(
            f"line {record.line_no}: exon strand ({record.strand}) != gene "
            f"strand ({gene_obj.strand}) for {gene_obj.id}"
        ),
    )
    exon = model_domain.Exon(
        record.start_pos,
        record.end_pos,
        record.chrom,
        strand,
        record.annots,
    )
    if gene_obj.add_exon(isoform, exon):
        ctx.stats.exon_count += 1
    if created_isoform:
        ctx.stats.iso_count += 1
    transcript = gene_obj.transcripts.get(isoform.id)
    if isinstance(transcript, model_domain.Transcript):
        ctx.model.mrna_forms.setdefault(record.chrom, {})[transcript.id] = transcript
        ctx.model.mrna_gene.setdefault(record.chrom, {})[transcript.id] = gene_obj
        _drain_pending_children(
            ctx,
            transcript=transcript,
            parent_gene=gene_obj,
            chrom=record.chrom,
        )


def handle_mrna_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not _has_chromosome(
        ctx,
        record,
        ctx.model.model,
        fail_message=(
            f"line {record.line_no}: Mrna with missing chromosome dictionary {record.chrom} "
            f"(known: {_known_chromosomes(ctx.model.model)})"
        ),
    ):
        return

    transcript_id = annotation_value(record.annots, model_domain.ID_FIELD)
    if transcript_id is None:
        ctx.fail(f"line {record.line_no}: Mrna with missing ID")
        return
    transcript_id = transcript_id.upper()

    parent_id = parent_annotation(record)
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

    strand = resolve_strand(
        ctx,
        observed=record.strand,
        expected=mrna_gene.strand,
        mismatch_message=(
            f"line {record.line_no}: Mrna strand ({record.strand}) does not "
            f"match gene strand ({mrna_gene.strand})"
        ),
    )
    mrna_attr = {
        str(model_domain.PARENT_FIELD): mrna_gene.id,
        str(model_domain.NAME_FIELD): transcript_id,
        str(model_domain.ID_FIELD): transcript_id,
    }
    mrna = model_domain.Transcript(
        transcript_id,
        record.start_pos,
        record.end_pos,
        record.chrom,
        strand,
        attr=mrna_attr,
    )
    mrna_gene.add_transcript(mrna)
    ctx.model.mrna_forms.setdefault(record.chrom, {})[transcript_id] = mrna
    ctx.model.mrna_gene.setdefault(record.chrom, {})[transcript_id] = mrna_gene
    _drain_pending_children(
        ctx,
        transcript=mrna,
        parent_gene=mrna_gene,
        chrom=record.chrom,
    )
    ctx.stats.mrna_count += 1


def handle_transcript_region_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not _has_chromosome(
        ctx,
        record,
        ctx.model.model,
        fail_message=(
            f"line {record.line_no}: {record.rec_type} has unrecognized chromosome: "
            f"{record.chrom} (known: {_known_chromosomes(ctx.model.model)})"
        ),
    ):
        return

    parent_id = parent_annotation(record)
    if parent_id is None:
        return
    mrna_record = resolve_parent(ctx, parent_id, record.chrom, search_genes=False)
    if mrna_record is None:
        ctx.write_verbose(f"line {record.line_no}: no Mrna {parent_id} found for {record.rec_type}")
        _queue_pending_child(ctx.pending_regions, record, parent_id)
        return
    if not isinstance(mrna_record, model_domain.Transcript):
        ctx.fail(f"line {record.line_no}: parent {parent_id} is not an Mrna record")
        return
    parent_gene = ctx.model.get_mrna_parent(record.chrom, mrna_record.id)
    if parent_gene is None:
        ctx.fail(
            f"line {record.line_no}: Mrna {mrna_record.id} is missing a parent gene for CDS record"
        )
        return

    strand = resolve_strand(
        ctx,
        observed=record.strand,
        expected=mrna_record.strand,
        mismatch_message=(
            f"line {record.line_no}: CDS strand ({record.strand}) does not "
            f"match Mrna strand ({mrna_record.strand})"
        ),
    )
    cds = model_domain.cds_factory(
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


def handle_misc_feature_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not _has_chromosome(ctx, record, ctx.model.model):
        return
    parent_value = parent_annotation(record)
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
            model_domain.BaseFeature(
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


def handle_chromosome_record(ctx: ParseContext, record: ParsedRecord) -> None:
    ctx.model.add_chromosome(record.start_pos, record.end_pos, record.chrom)


RECORD_HANDLERS: dict[RecordType, RecordHandler] = {
    RecordType.GENE: handle_gene_record,
    RecordType.PSEUDOGENE: handle_gene_record,
    RecordType.EXON: handle_exon_record,
    RecordType.PSEUDOGENIC_EXON: handle_exon_record,
    RecordType.MRNA: handle_mrna_record,
    RecordType.PSEUDOGENIC_TRANSCRIPT: handle_mrna_record,
    RecordType.CDS: handle_transcript_region_record,
    RecordType.FIVE_PRIME_UTR: handle_transcript_region_record,
    RecordType.THREE_PRIME_UTR: handle_transcript_region_record,
    RecordType.CHROMOSOME: handle_chromosome_record,
}


__all__ = [
    "RECORD_HANDLERS",
    "RecordHandler",
    "candidate_parent_ids",
    "handle_misc_feature_record",
    "resolve_parent",
]
