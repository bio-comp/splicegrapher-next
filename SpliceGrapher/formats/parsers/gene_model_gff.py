"""GFF parser boundary for ``GeneModel`` loading."""

from __future__ import annotations

from collections.abc import Sequence

import structlog

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.formats.parsers.gene_model_gff_context import GeneModelLike, ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_handlers import (
    RECORD_HANDLERS,
    handle_misc_feature_record,
)
from SpliceGrapher.formats.parsers.gene_model_gff_records import (
    iter_record_lines,
    normalize_chromosomes,
    parse_record_line,
)
from SpliceGrapher.shared.format_utils import comma_format
from SpliceGrapher.shared.progress import ProgressIndicator

LOGGER = structlog.get_logger(__name__)


def load_gene_model_records(
    model: GeneModelLike,
    gff_records: model_domain.GffRecordSource,
    *,
    require_notes: bool = False,
    chromosomes: Sequence[str] | str | None = None,
    verbose: bool = False,
    ignore_errors: bool = False,
) -> None:
    """Load gene model records from GFF-like input into ``model``."""
    chrom_filter = normalize_chromosomes(chromosomes)
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
    for line_no, line in enumerate(iter_record_lines(gff_records, verbose=verbose), start=1):
        indicator.update()
        record = parse_record_line(ctx, line, line_no)
        if record is None:
            continue
        handler = RECORD_HANDLERS.get(record.rec_type, handle_misc_feature_record)
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
                parent_cache_clear_count=ctx.stats.parent_cache_clear_count,
                parent_candidate_clear_count=ctx.stats.parent_candidate_clear_count,
            )
        else:
            LOGGER.warning("gene_model_gff_no_genes_loaded")


__all__ = ["load_gene_model_records"]
