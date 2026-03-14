"""Parent, strand, and isoform resolution helpers for gene-model GFF parsing."""

from __future__ import annotations

from collections.abc import Container, Mapping

import SpliceGrapher.formats.models as model_domain
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_records import (
    ParsedRecord,
    annotation_value,
)


def _subnames(full_string: str, delimiter: str) -> list[str]:
    parts = full_string.split(delimiter)
    return [delimiter.join(parts[:i]) for i in range(len(parts) - 1, 0, -1)]


def known_chromosomes(mapping: Mapping[str, Mapping[str, model_domain.Gene]]) -> str:
    return ",".join(mapping.keys())


def has_chromosome(
    ctx: ParseContext,
    record: ParsedRecord,
    mapping: Container[str],
    *,
    fail_message: str | None = None,
) -> bool:
    if record.chrom in mapping:
        return True
    if fail_message is not None:
        ctx.fail(fail_message)
    return False


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


def resolve_exon_parent(ctx: ParseContext, record: ParsedRecord) -> model_domain.Gene | None:
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


def resolve_isoform(
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


def resolve_mrna_gene(
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


def extract_parent_gene_id(parent_id: str) -> str | None:
    if not parent_id:
        return None
    primary_parent = parent_id.split(",", 1)[0].strip()
    if not primary_parent:
        return None
    return primary_parent.split(".", 1)[0]


__all__ = [
    "candidate_parent_ids",
    "extract_parent_gene_id",
    "has_chromosome",
    "known_chromosomes",
    "resolve_exon_parent",
    "resolve_isoform",
    "resolve_mrna_gene",
    "resolve_parent",
    "resolve_strand",
]
