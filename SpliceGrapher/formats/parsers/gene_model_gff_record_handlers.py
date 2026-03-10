"""Record dispatch handlers for gene-model GFF loading."""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeAlias

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_pending import (
    drain_pending_children,
    queue_pending_child,
)
from SpliceGrapher.formats.parsers.gene_model_gff_records import (
    ParsedRecord,
    annotation_value,
    clean_name,
    parent_annotation,
)
from SpliceGrapher.formats.parsers.gene_model_gff_resolution import (
    extract_parent_gene_id,
    has_chromosome,
    known_chromosomes,
    resolve_exon_parent,
    resolve_isoform,
    resolve_mrna_gene,
    resolve_parent,
    resolve_strand,
)

RecordHandler: TypeAlias = Callable[[ParseContext, ParsedRecord], None]


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
    if not has_chromosome(ctx, record, ctx.model.model):
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

    gene_obj = resolve_exon_parent(ctx, record)
    if gene_obj is None:
        queue_pending_child(ctx.pending_exons, record, parent_id)
        return

    isoform, created_isoform = resolve_isoform(ctx, record, gene_obj)
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
        drain_pending_children(
            ctx,
            transcript=transcript,
            parent_gene=gene_obj,
            chrom=record.chrom,
        )


def handle_mrna_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not has_chromosome(
        ctx,
        record,
        ctx.model.model,
        fail_message=(
            f"line {record.line_no}: Mrna with missing chromosome dictionary {record.chrom} "
            f"(known: {known_chromosomes(ctx.model.model)})"
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

    mrna_gene = resolve_mrna_gene(ctx, record, parent_id)
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
    drain_pending_children(
        ctx,
        transcript=mrna,
        parent_gene=mrna_gene,
        chrom=record.chrom,
    )
    ctx.stats.mrna_count += 1


def handle_transcript_region_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not has_chromosome(
        ctx,
        record,
        ctx.model.model,
        fail_message=(
            f"line {record.line_no}: {record.rec_type} has unrecognized chromosome: "
            f"{record.chrom} (known: {known_chromosomes(ctx.model.model)})"
        ),
    ):
        return

    parent_id = parent_annotation(record)
    if parent_id is None:
        return
    mrna_record = resolve_parent(ctx, parent_id, record.chrom, search_genes=False)
    if mrna_record is None:
        ctx.write_verbose(f"line {record.line_no}: no Mrna {parent_id} found for {record.rec_type}")
        queue_pending_child(ctx.pending_regions, record, parent_id)
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


def handle_misc_feature_record(ctx: ParseContext, record: ParsedRecord) -> None:
    if not has_chromosome(ctx, record, ctx.model.model):
        return
    parent_value = parent_annotation(record)
    if parent_value is None:
        return
    parent_gene_id = extract_parent_gene_id(parent_value)
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


__all__ = ["RECORD_HANDLERS", "RecordHandler", "handle_misc_feature_record"]
