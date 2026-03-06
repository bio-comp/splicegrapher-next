from __future__ import annotations

from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats import gene_model as model_facade
from SpliceGrapher.formats import models as model_domain
from SpliceGrapher.formats import serializers as ser
from SpliceGrapher.formats.gene_model import CDS, Exon, Gene, Transcript


def test_models_module_no_longer_owns_dict_text_serializers() -> None:
    assert not hasattr(model_domain, "dict_to_gff")
    assert not hasattr(model_domain, "dict_to_gtf")


def test_gene_model_facade_no_longer_reexports_dict_text_serializers() -> None:
    assert not hasattr(model_facade, "dict_to_gff")
    assert not hasattr(model_facade, "dict_to_gtf")


def test_models_and_facade_no_longer_export_mrna_isoform_aliases() -> None:
    assert not hasattr(model_domain, "Mrna")
    assert not hasattr(model_domain, "Isoform")
    assert not hasattr(model_facade, "Mrna")
    assert not hasattr(model_facade, "Isoform")


def test_transcript_gtf_lines_remain_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    transcript = Transcript("TX1", 1, 200, "chr1", "-", feature_type=RecordType.MRNA)
    transcript.parent = gene.id
    transcript.add_cds(CDS(80, 90, "chr1", "-"))
    transcript.add_cds(CDS(20, 30, "chr1", "-"))

    cds_starts = [
        int(line.split("\t")[3])
        for line in ser.gtf_lines_for_transcript(transcript, gene)
        if line.split("\t")[2] == RecordType.CDS.value
    ]
    assert cds_starts == [20, 80]


def test_gene_gtf_text_common_path_remains_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    isoform = Transcript("TX1", 1, 200, "chr1", "-")
    transcript = Transcript("TX1", 1, 200, "chr1", "-", feature_type=RecordType.MRNA)
    gene.add_transcript(isoform)
    gene.add_transcript(transcript)

    isoform.add_exon(Exon(70, 75, "chr1", "-"))
    isoform.add_exon(Exon(20, 30, "chr1", "-"))
    transcript.add_cds(CDS(80, 90, "chr1", "-"))
    transcript.add_cds(CDS(25, 35, "chr1", "-"))

    feature_starts = [
        int(line.split("\t")[3])
        for line in ser.gtf_text_for_gene(gene).splitlines()
        if line.split("\t")[2] in {RecordType.EXON.value, RecordType.CDS.value}
    ]
    assert feature_starts == sorted(feature_starts)
