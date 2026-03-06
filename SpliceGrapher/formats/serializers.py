"""Pure formatting helpers for gene-model domain entities.

This module intentionally owns text serialization concerns (GFF/GTF) so the
domain model can stay focused on biological state and query behavior.
"""

from __future__ import annotations

from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats import models as model_domain


def dict_to_gff(d: dict[str, str]) -> str:
    """Render an attribute mapping using GFF3 key=value pairs."""
    return ";".join([f"{k}={v}" for k, v in sorted(d.items()) if k != "parent"])


def dict_to_gtf(d: dict[str, str]) -> str:
    """Render an attribute mapping using GTF key \"value\" pairs."""
    return "; ".join([f'{k} "{v}"' for k, v in sorted(d.items()) if k != "parent"])


def gff_line_for_chromosome(chromosome: model_domain.Chromosome) -> str:
    name_str = chromosome.name.capitalize()
    return (
        f"{name_str}\t{model_domain.GFF_ID}\tchromosome\t"
        f"{chromosome.start()}\t{chromosome.end()}\t.\t.\t.\tID={name_str};Name={name_str}"
    )


def gff_line_for_feature(
    feature: model_domain.BaseFeature | model_domain.Transcript | model_domain.Gene,
    alt_attributes: model_domain.ExtraAttrMap | None = None,
) -> str:
    attrs = dict(feature.attributes)
    if alt_attributes is not None:
        attrs.update({str(k): str(v) for k, v in alt_attributes.items()})
    return (
        f"{feature.chromosome.capitalize()}\t{model_domain.GFF_ID}\t{feature.feature_type}\t"
        f"{feature.minpos}\t{feature.maxpos}\t.\t{feature.strand}\t.\t"
        f"{dict_to_gff(attrs)}"
    )


def gtf_line_for_feature(
    feature: model_domain.BaseFeature,
    transcript_id: str,
    gene: model_domain.Gene,
    exon_id: int,
) -> str:
    attrs = {str(k).lower(): str(v) for k, v in feature.attributes.items()}
    attrs[model_domain.GTF_GENE_ID] = gene.id
    attrs[model_domain.GTF_GENE_NAME] = gene.name if gene.name is not None else gene.id
    attrs[model_domain.GTF_TRANSCRIPT] = transcript_id
    attrs[model_domain.GTF_EXON_ID] = str(exon_id)
    source = attrs.get(model_domain.GTF_SOURCE, model_domain.GFF_ID)
    return (
        f"{feature.chromosome.capitalize()}\t{source}\t{feature.feature_type}\t"
        f"{feature.minpos}\t{feature.maxpos}\t.\t{feature.strand}\t.\t"
        f"{dict_to_gtf(attrs)}"
    )


def gff_lines_for_transcript(transcript: model_domain.Transcript) -> list[str]:
    result = [gff_line_for_feature(transcript)]
    gff_attr = {model_domain.PARENT_FIELD: transcript.id}
    gtf_attr = {
        model_domain.PARENT_FIELD: transcript.id,
        model_domain.GTF_TRANSCRIPT: transcript.id,
        model_domain.GTF_TRANSNAME: transcript.id,
        model_domain.GTF_PROTEIN_ID: transcript.id,
    }
    all_children: list[model_domain.BaseFeature] = list(transcript.exons) + list(transcript.cds)
    all_children.sort(key=model_domain.feature_sort_key)
    for child in all_children:
        attrs = gtf_attr if model_domain.GTF_TRANSCRIPT in child.attributes else gff_attr
        result.append(gff_line_for_feature(child, alt_attributes=attrs))
    return result


def gtf_start_codon_line(
    transcript: model_domain.Transcript,
    gene: model_domain.Gene,
) -> str | None:
    if transcript.start_codon_pos is None:
        return None
    return (
        f"{transcript.chromosome}\t{model_domain.GFF_ID}\tstart_codon\t"
        f"{transcript.start_codon_pos[0]}\t{transcript.start_codon_pos[1]}\t.\t"
        f'{transcript.strand}\t.\tgene_id "{gene.id}"; transcript_id "{transcript.id}"'
    )


def gtf_stop_codon_line(
    transcript: model_domain.Transcript,
    gene: model_domain.Gene,
) -> str | None:
    if transcript.end_codon_pos is None:
        return None
    return (
        f"{transcript.chromosome}\t{model_domain.GFF_ID}\tstop_codon\t"
        f"{transcript.end_codon_pos[0]}\t{transcript.end_codon_pos[1]}\t.\t"
        f'{transcript.strand}\t.\tgene_id "{gene.id}"; transcript_id "{transcript.id}"'
    )


def gtf_lines_for_transcript(
    transcript: model_domain.Transcript,
    gene: model_domain.Gene,
) -> list[str]:
    result: list[str] = []
    first_codon = (
        gtf_start_codon_line(transcript, gene)
        if transcript.strand == "+"
        else gtf_stop_codon_line(transcript, gene)
    )
    if first_codon:
        result.append(first_codon)

    all_features: list[model_domain.BaseFeature] = list(transcript.exons) + list(transcript.cds)
    all_features.sort(key=model_domain.gtf_feature_sort_key)
    exon_counter = 0
    cds_counter = 0
    for feature in all_features:
        if feature.feature_type == RecordType.EXON:
            exon_counter += 1
            result.append(gtf_line_for_feature(feature, transcript.id, gene, exon_counter))
        elif feature.feature_type == RecordType.CDS:
            cds_counter += 1
            result.append(gtf_line_for_feature(feature, transcript.id, gene, cds_counter))

    second_codon = (
        gtf_stop_codon_line(transcript, gene)
        if transcript.strand == "+"
        else gtf_start_codon_line(transcript, gene)
    )
    if second_codon:
        result.append(second_codon)
    return result


def gff_text_for_gene(gene: model_domain.Gene) -> str:
    lines = [gff_line_for_feature(gene, alt_attributes={model_domain.ID_FIELD: gene.id})]
    for transcript_id in sorted(gene.transcripts):
        lines.extend(gff_lines_for_transcript(gene.transcripts[transcript_id]))
    return "\n".join(lines)


def gtf_text_for_gene(gene: model_domain.Gene) -> str:
    lines: list[str] = []
    for transcript_id in sorted(gene.transcripts):
        lines.extend(gtf_lines_for_transcript(gene.transcripts[transcript_id], gene))
    return "\n".join(lines)


__all__ = [
    "dict_to_gff",
    "dict_to_gtf",
    "gff_line_for_chromosome",
    "gff_line_for_feature",
    "gff_lines_for_transcript",
    "gff_text_for_gene",
    "gtf_line_for_feature",
    "gtf_lines_for_transcript",
    "gtf_start_codon_line",
    "gtf_stop_codon_line",
    "gtf_text_for_gene",
]
