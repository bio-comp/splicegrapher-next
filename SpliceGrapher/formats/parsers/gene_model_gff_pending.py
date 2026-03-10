"""Pending child queue helpers for gene-model GFF parsing."""

from __future__ import annotations

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_records import ParsedRecord

from .gene_model_gff_resolution import resolve_strand


def pending_key(chrom: str, transcript_id: str) -> tuple[str, str]:
    return (chrom, transcript_id.upper())


def queue_pending_child(
    pending_map: dict[tuple[str, str], list[ParsedRecord]],
    record: ParsedRecord,
    parent_id: str,
) -> None:
    pending_map.setdefault(pending_key(record.chrom, parent_id), []).append(record)


def drain_pending_children(
    ctx: ParseContext,
    *,
    transcript: model_domain.Transcript,
    parent_gene: model_domain.Gene,
    chrom: str,
) -> None:
    key = pending_key(chrom, transcript.id)

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


__all__ = ["drain_pending_children", "pending_key", "queue_pending_child"]
