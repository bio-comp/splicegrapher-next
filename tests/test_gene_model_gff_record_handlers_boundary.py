from __future__ import annotations

import importlib

from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats import gene_model as gm
from SpliceGrapher.formats.parsers.gene_model_gff_context import ParseContext
from SpliceGrapher.formats.parsers.gene_model_gff_records import ParsedRecord

parser_record_handlers = importlib.import_module(
    "SpliceGrapher.formats.parsers.gene_model_gff_record_handlers"
)


def _build_context() -> tuple[ParseContext, gm.Gene]:
    model = gm.GeneModel()
    model.add_chromosome(1, 500, "chr1")
    gene = gm.Gene("GENE1", None, 10, 200, "chr1", "+")
    model.add_gene(gene)
    ctx = ParseContext(
        model=model,
        require_notes=False,
        verbose=False,
        ignore_errors=False,
    )
    return ctx, gene


def _build_record(*, parent_id: str) -> ParsedRecord:
    return ParsedRecord(
        line_no=1,
        raw_line=f"chr1\tsrc\tcds\t10\t20\t.\t+\t.\tParent={parent_id}",
        annots={str(gm.PARENT_FIELD): parent_id},
        chrom="chr1",
        rec_type=RecordType.CDS,
        start_pos=10,
        end_pos=20,
        strand="+",
    )


def _build_transcript(transcript_id: str, parent_gene_id: str) -> gm.Transcript:
    return gm.Transcript(
        transcript_id,
        10,
        200,
        "chr1",
        "+",
        attr={
            str(gm.PARENT_FIELD): parent_gene_id,
            str(gm.ID_FIELD): transcript_id,
            str(gm.NAME_FIELD): transcript_id,
        },
    )


def test_register_transcript_links_tracks_mrna_maps() -> None:
    ctx, gene = _build_context()
    transcript = _build_transcript("TX1", gene.id)

    parser_record_handlers._register_transcript_links(
        ctx,
        transcript=transcript,
        parent_gene=gene,
        chrom="chr1",
    )

    assert ctx.model.mrna_forms["chr1"]["TX1"] is transcript
    assert ctx.model.mrna_gene["chr1"]["TX1"] is gene


def test_resolve_region_target_returns_transcript_and_parent_gene_when_present() -> None:
    ctx, gene = _build_context()
    transcript = _build_transcript("TX1", gene.id)
    gene.add_transcript(transcript)
    ctx.model.mrna_forms.setdefault("chr1", {})["TX1"] = transcript
    ctx.model.mrna_gene.setdefault("chr1", {})["TX1"] = gene
    record = _build_record(parent_id="TX1")

    result = parser_record_handlers._resolve_region_target(ctx, record, "TX1")

    assert result is not None
    resolved_transcript, parent_gene = result
    assert resolved_transcript is transcript
    assert parent_gene is gene


def test_resolve_region_target_queues_pending_region_when_parent_missing() -> None:
    ctx, _ = _build_context()
    record = _build_record(parent_id="MISSING_TX")

    result = parser_record_handlers._resolve_region_target(ctx, record, "MISSING_TX")

    assert result is None
    assert ("chr1", "MISSING_TX") in ctx.pending_regions
