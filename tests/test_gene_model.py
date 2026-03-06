from __future__ import annotations

import inspect
from collections.abc import Sequence
from pathlib import Path

import pytest

from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.core.interval_helpers import InMemoryIntervalIndex
from SpliceGrapher.formats import gene_model as gm
from SpliceGrapher.formats import serializers as ser
from SpliceGrapher.formats.gene_model import Exon, Gene, GeneModel


def test_interval_index_search_supports_recursive_midpoint_indexing() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert (
        InMemoryIntervalIndex(features).predecessor_or_containing(query) is features[1]
    )


def test_mrna_sorted_exons_uses_explicit_minintron_keyword() -> None:
    transcript = gm.Transcript("tx1", 1, 20, "chr1", "+", feature_type=RecordType.MRNA)
    transcript.add_cds(gm.CDS(1, 5, "chr1", "+"))
    transcript.add_cds(gm.CDS(7, 10, "chr1", "+"))

    merged = transcript.sorted_exons(minintron=3)
    unmerged = transcript.sorted_exons(minintron=2)

    assert len(merged) == 1
    assert len(unmerged) == 2

    with pytest.raises(TypeError):
        transcript.sorted_exons(min_intron=3)  # type: ignore[call-arg]


def test_interval_index_search_returns_expected_match() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert (
        InMemoryIntervalIndex(features).predecessor_or_containing(query) is features[1]
    )


def test_feature_overlap_and_contains_match_interval_helper_semantics() -> None:
    left = Exon(10, 20, "chr1", "+")
    right = Exon(20, 30, "chr1", "+")
    nested = Exon(12, 18, "chr1", "+")

    assert gm.feature_overlaps(left, right)
    assert gm.feature_contains(left, nested)
    assert not gm.feature_contains(left, right)


def test_interval_index_search_returns_preceding_feature_when_not_contained() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(32, 33, "chr1", "+")

    assert (
        InMemoryIntervalIndex(features).predecessor_or_containing(query) is features[1]
    )


def test_interval_index_search_rejects_invalid_bounds() -> None:
    features = [Exon(1, 10, "chr1", "+"), Exon(20, 30, "chr1", "+")]
    query = Exon(8, 9, "chr1", "+")

    with pytest.raises(ValueError, match="Invalid search bounds"):
        InMemoryIntervalIndex(features).predecessor_or_containing(query, lo=2, hi=1)


def test_get_gene_from_locations_uses_interval_index() -> None:
    model = GeneModel()
    model.add_chromosome(1, 1000, "chr1")
    gene = Gene("GENE_A", None, 150, 190, "chr1", "+", name="GENE_A")
    model.add_gene(gene)
    model.make_sorted_model()

    assert model.get_gene_from_locations("chr1", 155, 160, "+") is gene


def test_get_genes_in_range_uses_interval_index_when_available(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = GeneModel()
    model.add_chromosome(1, 1000, "chr1")
    gene = Gene("GENE_A", None, 150, 190, "chr1", "+", name="GENE_A")
    model.add_gene(gene)
    model.make_sorted_model()

    def _fail_scan(*args: object, **kwargs: object) -> list[Gene]:
        raise AssertionError(
            "get_gene_records scan path should not be used when index exists"
        )

    monkeypatch.setattr(GeneModel, "get_gene_records", _fail_scan)

    hits = model.get_genes_in_range("chr1", 100, 200)

    assert hits == [gene]


def test_gene_model_init_rejects_implicit_autoload() -> None:
    with pytest.raises(TypeError, match="positional"):
        GeneModel("fake.gff3")  # type: ignore[misc,arg-type]


def test_gene_model_constructor_no_longer_exposes_gff_path_parameter() -> None:
    signature = inspect.signature(GeneModel)
    assert "gff_path" not in signature.parameters


def test_gene_model_from_gff_factory_loads_records() -> None:
    records = ["chr1\tsource\tgene\t10\t20\t.\t+\t.\tID=GENE1;Name=GENE1"]
    model = GeneModel.from_gff(records)

    assert model.get_gene("chr1", "GENE1") is not None


def test_parser_collapses_exon_and_cds_into_single_transcript_object() -> None:
    records = [
        "chr1\tsource\tgene\t1\t200\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsource\tmrna\t1\t200\t.\t+\t.\tID=TX1;Parent=GENE1",
        "chr1\tsource\texon\t1\t50\t.\t+\t.\tID=EX1;Parent=TX1",
        "chr1\tsource\tcds\t10\t40\t.\t+\t0\tID=CDS1;Parent=TX1",
    ]
    model = GeneModel.from_gff(records)
    gene = model.get_gene("chr1", "GENE1")

    assert gene is not None
    assert "TX1" in gene.transcripts
    transcript = gene.transcripts["TX1"]
    assert len(transcript.exons) == 1
    assert len(transcript.cds) == 1
    assert gene.transcripts["TX1"] is transcript


def test_parser_unified_transcript_clean_order_is_singular_complete_and_codon_ready() -> (
    None
):
    records = [
        "chr1\tsource\tgene\t1\t200\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsource\tmrna\t1\t200\t.\t+\t.\tID=TX1;Parent=GENE1",
        "chr1\tsource\texon\t1\t30\t.\t+\t.\tID=EX1;Parent=TX1",
        "chr1\tsource\texon\t40\t70\t.\t+\t.\tID=EX2;Parent=TX1",
        "chr1\tsource\tcds\t10\t20\t.\t+\t0\tID=CDS1;Parent=TX1",
        "chr1\tsource\tcds\t50\t60\t.\t+\t0\tID=CDS2;Parent=TX1",
    ]
    model = GeneModel.from_gff(records)
    gene = model.get_gene("chr1", "GENE1")

    assert gene is not None
    assert len(gene.transcripts) == 1
    transcript = gene.transcripts["TX1"]
    assert len(transcript.exons) > 0
    assert len(transcript.cds) > 0

    transcript.infer_codons()
    assert transcript.start_codon_pos == (10, 12)
    assert transcript.end_codon_pos == (58, 60)


def test_parser_unified_transcript_handles_out_of_order_exon_and_cds() -> None:
    records = [
        "chr1\tsource\tgene\t1\t200\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsource\texon\t1\t50\t.\t+\t.\tID=EX1;Parent=TX1",
        "chr1\tsource\tcds\t10\t40\t.\t+\t0\tID=CDS1;Parent=TX1",
        "chr1\tsource\tmrna\t1\t200\t.\t+\t.\tID=TX1;Parent=GENE1",
    ]
    model = GeneModel.from_gff(records)
    gene = model.get_gene("chr1", "GENE1")

    assert gene is not None
    assert len(gene.transcripts) == 1
    transcript = gene.transcripts["TX1"]
    assert len(transcript.exons) == 1
    assert len(transcript.cds) == 1


def test_explicit_exons_drive_intron_and_junction_counts() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    gene.add_transcript(transcript)
    transcript.add_exon(Exon(1, 10, "chr1", "+"))
    transcript.add_exon(Exon(21, 30, "chr1", "+"))
    transcript.add_exon(Exon(41, 50, "chr1", "+"))
    transcript.add_cds(gm.CDS(2, 9, "chr1", "+"))
    transcript.add_cds(gm.CDS(22, 29, "chr1", "+"))
    transcript.add_cds(gm.CDS(42, 49, "chr1", "+"))

    assert len(gene.get_introns()) == 2
    assert len(gene.get_junctions()) == 2
    assert len(gene.acceptor_list()) == 2
    assert len(gene.donor_list()) == 2


def test_get_junctions_uses_inferred_exons_when_transcript_has_only_cds() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    gene.add_transcript(transcript)
    transcript.add_cds(gm.CDS(1, 10, "chr1", "+"))
    transcript.add_cds(gm.CDS(21, 30, "chr1", "+"))
    transcript.add_cds(gm.CDS(41, 50, "chr1", "+"))

    assert gene.get_introns() == {(10, 21), (30, 41)}
    assert gene.get_junctions() == {(10, 19), (30, 39)}


def test_transcript_sorted_exons_prefers_explicit_exons_over_cds_inference() -> None:
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    transcript.add_exon(Exon(1, 10, "chr1", "+"))
    transcript.add_exon(Exon(21, 30, "chr1", "+"))
    transcript.add_cds(gm.CDS(1, 30, "chr1", "+"))

    exons = transcript.sorted_exons()
    assert [(exon.minpos, exon.maxpos) for exon in exons] == [(1, 10), (21, 30)]


def test_transcript_sorted_exons_does_not_merge_overlapping_cds_regions() -> None:
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    transcript.add_cds(gm.CDS(10, 20, "chr1", "+"))
    transcript.add_cds(gm.CDS(19, 30, "chr1", "+"))

    exons = transcript.sorted_exons(minintron=2)

    assert [(exon.minpos, exon.maxpos) for exon in exons] == [(10, 20), (19, 30)]


def test_gene_is_single_exon_requires_all_transcripts_to_be_single_exon() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    tx1 = gm.Transcript("TX1", 1, 100, "chr1", "+")
    tx2 = gm.Transcript("TX2", 1, 100, "chr1", "+")

    gene.add_exon(tx1, Exon(1, 10, "chr1", "+"))
    gene.add_exon(tx2, Exon(20, 30, "chr1", "+"))
    assert gene.is_single_exon()

    gene.add_exon(tx2, Exon(40, 50, "chr1", "+"))
    assert not gene.is_single_exon()


def test_add_exon_keeps_existing_record_and_merges_new_annotation_keys() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    first = Exon(10, 20, "chr1", "+", {gm.ID_FIELD: "EX1", "source": "manual"})
    second = Exon(10, 20, "chr1", "+", {gm.ID_FIELD: "EX2", "evidence": "rna_seq"})

    added_first = gene.add_exon(transcript, first)
    added_second = gene.add_exon(transcript, second)

    stored = gene.exon_map[(10, 20)]
    assert added_first
    assert not added_second
    assert stored.attributes[str(gm.ID_FIELD)] == "EX1"
    assert stored.attributes["source"] == "manual"
    assert stored.attributes["evidence"] == "rna_seq"


def test_add_cds_keeps_existing_record_and_merges_new_annotation_keys() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    first = gm.CDS(10, 20, "chr1", "+", {gm.ID_FIELD: "CDS1", "phase": "0"})
    second = gm.CDS(10, 20, "chr1", "+", {gm.ID_FIELD: "CDS2", "source": "rna_seq"})

    added_first = gene.add_cds(transcript, first)
    added_second = gene.add_cds(transcript, second)

    stored = gene.cds_map[(RecordType.CDS, 10, 20)]
    assert added_first
    assert not added_second
    assert stored.attributes[str(gm.ID_FIELD)] == "CDS1"
    assert stored.attributes["phase"] == "0"
    assert stored.attributes["source"] == "rna_seq"


def test_transcript_add_exon_merges_non_conflicting_attributes_on_duplicates() -> None:
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    first = Exon(10, 20, "chr1", "+", {"source": "manual"})
    second = Exon(10, 20, "chr1", "+", {"evidence": "rna_seq"})

    assert transcript.add_exon(first)
    assert not transcript.add_exon(second)

    stored = transcript.exon_map[(10, 20)]
    assert stored.attributes["source"] == "manual"
    assert stored.attributes["evidence"] == "rna_seq"


def test_transcript_add_cds_merges_non_conflicting_attributes_on_duplicates() -> None:
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    first = gm.CDS(10, 20, "chr1", "+", {"phase": "0"})
    second = gm.CDS(10, 20, "chr1", "+", {"source": "rna_seq"})

    assert transcript.add_cds(first)
    assert not transcript.add_cds(second)

    stored = transcript.cds_map[(RecordType.CDS, 10, 20)]
    assert stored.attributes["phase"] == "0"
    assert stored.attributes["source"] == "rna_seq"


def test_transcript_add_cds_treats_feature_type_as_part_of_dedupe_key() -> None:
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    coding = gm.CDS(10, 20, "chr1", "+", {"id": "coding"})
    utr = gm.FpUtr(10, 20, "chr1", "+", {"id": "utr"})

    assert transcript.add_cds(coding)
    assert transcript.add_cds(utr)
    assert len(transcript.cds) == 2


def test_add_transcript_registers_preloaded_exons_and_cds_in_gene_catalogs() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 100, "chr1", "+")
    transcript.add_exon(Exon(10, 20, "chr1", "+", {"source": "manual"}))
    transcript.add_cds(gm.CDS(12, 18, "chr1", "+", {"phase": "0"}))

    gene.add_transcript(transcript)

    assert len(gene.exons) == 1
    assert len(gene.cds) == 1
    assert gene.exon_map[(10, 20)].attributes["source"] == "manual"
    assert gene.cds_map[(RecordType.CDS, 12, 18)].attributes["phase"] == "0"


def test_add_transcript_dedupes_gene_exons_across_transcript_ids() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    tx1 = gm.Transcript("TX1", 1, 100, "chr1", "+")
    tx2 = gm.Transcript("TX2", 1, 100, "chr1", "+")
    tx1.add_exon(Exon(10, 20, "chr1", "+", {"source": "manual"}))
    tx2.add_exon(Exon(10, 20, "chr1", "+", {"evidence": "rna_seq"}))

    gene.add_transcript(tx1)
    gene.add_transcript(tx2)

    assert len(gene.exons) == 1
    merged = gene.exon_map[(10, 20)]
    assert merged.attributes["source"] == "manual"
    assert merged.attributes["evidence"] == "rna_seq"


def test_add_transcript_rejects_chromosome_mismatch_without_mutation() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    bad_tx = gm.Transcript("TX1", 1, 100, "chr2", "+")

    with pytest.raises(ValueError, match="chromosome"):
        gene.add_transcript(bad_tx)

    assert not gene.transcripts
    assert not gene.exons
    assert not gene.cds


def test_add_transcript_rejects_strand_mismatch_without_mutation() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    bad_tx = gm.Transcript("TX1", 1, 100, "chr1", "-")

    with pytest.raises(ValueError, match="strand"):
        gene.add_transcript(bad_tx)

    assert not gene.transcripts
    assert not gene.exons
    assert not gene.cds


def test_add_exon_mismatched_transcript_does_not_partially_mutate_gene() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    bad_tx = gm.Transcript("TX1", 1, 100, "chr2", "+")
    exon = Exon(10, 20, "chr1", "+")

    with pytest.raises(ValueError, match="chromosome"):
        gene.add_exon(bad_tx, exon)

    assert not gene.transcripts
    assert not gene.exons
    assert not gene.exon_map


def test_add_cds_mismatched_transcript_does_not_partially_mutate_gene() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+", name="GENE1")
    bad_tx = gm.Transcript("TX1", 1, 100, "chr2", "+")
    cds = gm.CDS(10, 20, "chr1", "+")

    with pytest.raises(ValueError, match="chromosome"):
        gene.add_cds(bad_tx, cds)

    assert not gene.transcripts
    assert not gene.cds
    assert not gene.cds_map


def test_write_gff_writes_chromosome_records(tmp_path: Path) -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    model.add_gene(
        Gene(
            id="GENE1",
            note=None,
            start_pos=10,
            end_pos=20,
            chromosome="chr1",
            strand="+",
            name="GENE1",
            attr={},
        )
    )

    gff_path = tmp_path / "model.gff3"
    model.write_gff(str(gff_path))

    lines = gff_path.read_text(encoding="utf-8").splitlines()
    expected_prefix = (
        "Chr1\tSpliceGrapher\tchromosome\t1\t100\t.\t.\t.\tID=Chr1;Name=Chr1"
    )
    assert lines[0].startswith(expected_prefix)


def test_get_annotation_dict_ignores_malformed_tokens_and_splits_once() -> None:
    model = GeneModel()
    parsed = model.get_annotation_dict("ID=GENE1;badtoken;Name=A=B;Parent=P1;=junk;")

    assert parsed == {"ID": "GENE1", "Name": "A=B", "Parent": "P1"}


def test_load_gene_model_rejects_unknown_record_type() -> None:
    model = GeneModel()
    records = ["chr1\tsource\tunknown_type\t1\t10\t.\t+\t.\tID=U1;Name=U1"]

    with pytest.raises(ValueError):
        model.load_gene_model(records)


def test_gene_model_record_type_collections_are_enum_backed() -> None:
    assert all(isinstance(item, RecordType) for item in gm.KNOWN_RECTYPES)
    assert all(isinstance(item, RecordType) for item in gm.IGNORE_RECTYPES)
    assert all(isinstance(item, RecordType) for item in gm.CDS_TYPES)


def test_gene_model_alias_mapping_uses_recordtype_domain() -> None:
    assert gm.RECTYPE_MAP[RecordType.PREDICTED_GENE] is RecordType.GENE
    assert gm.RECTYPE_MAP[RecordType.CDS_PREDICTED] is RecordType.CDS


def test_gene_model_valid_strands_uses_enum_domain() -> None:
    assert set(Strand) == gm.VALID_STRANDS


def test_gene_model_gtf_policy_and_splice_offset_contract_constants() -> None:
    assert gm.GTF_ORDER_POLICY == "genomic_ascending"
    assert gm.SPLICE_DIMER_OFFSET == 2


def test_exon_splice_site_offset_contract_is_explicit_for_both_strands() -> None:
    plus_exon = Exon(10, 20, "chr1", "+")
    minus_exon = Exon(10, 20, "chr1", "-")

    assert plus_exon.acceptor() == 10 - gm.SPLICE_DIMER_OFFSET
    assert plus_exon.donor() == 20
    assert minus_exon.acceptor() == 20
    assert minus_exon.donor() == 10 - gm.SPLICE_DIMER_OFFSET


def test_isoform_gtf_strings_remain_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    isoform = gm.Transcript("TX1", 1, 200, "chr1", "-", feature_type=gm.ISOFORM_TYPE)
    isoform.set_parent(gene.id)
    isoform.add_exon(Exon(80, 90, "chr1", "-"))
    isoform.add_exon(Exon(20, 30, "chr1", "-"))

    starts = [
        int(line.split("\t")[3]) for line in ser.gtf_lines_for_transcript(isoform, gene)
    ]
    assert starts == [20, 80]


def test_mrna_gtf_strings_remain_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    transcript = gm.Transcript("TX1", 1, 200, "chr1", "-", feature_type=RecordType.MRNA)
    transcript.parent = gene.id
    transcript.add_cds(gm.CDS(80, 90, "chr1", "-"))
    transcript.add_cds(gm.CDS(20, 30, "chr1", "-"))

    cds_starts = [
        int(line.split("\t")[3])
        for line in ser.gtf_lines_for_transcript(transcript, gene)
    ]
    assert cds_starts == [20, 80]


def test_gene_gtf_strings_common_path_remains_genomic_ascending_on_minus_strand() -> (
    None
):
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    isoform = gm.Transcript("TX1", 1, 200, "chr1", "-", feature_type=gm.ISOFORM_TYPE)
    transcript = gm.Transcript("TX1", 1, 200, "chr1", "-", feature_type=RecordType.MRNA)
    gene.add_transcript(isoform)
    gene.add_transcript(transcript)

    isoform.add_exon(Exon(70, 75, "chr1", "-"))
    isoform.add_exon(Exon(20, 30, "chr1", "-"))
    transcript.add_cds(gm.CDS(80, 90, "chr1", "-"))
    transcript.add_cds(gm.CDS(25, 35, "chr1", "-"))

    feature_starts = [
        int(line.split("\t")[3])
        for line in ser.gtf_text_for_gene(gene).splitlines()
        if line.split("\t")[2] in {RecordType.EXON.value, RecordType.CDS.value}
    ]
    assert feature_starts == sorted(feature_starts)


def test_exon_default_attributes_are_not_shared() -> None:
    exon_one = Exon(1, 10, "chr1", "+")
    exon_two = Exon(20, 30, "chr1", "+")

    exon_one.attributes["tag"] = "one"

    assert "tag" not in exon_two.attributes


def test_transcript_regions_do_not_store_parent_backlinks() -> None:
    gene = Gene("GENE1", None, 1, 100, "chr1", "+")
    isoform = gm.Transcript("TX1", 1, 100, "chr1", "+", feature_type=gm.ISOFORM_TYPE)
    exon = Exon(10, 20, "chr1", "+")

    gene.add_transcript(isoform)
    isoform.add_exon(exon)

    assert not hasattr(exon, "parents")


def test_model_classes_define_slots_without_instance_dict() -> None:
    isoform = gm.Transcript("TX1", 1, 20, "chr1", "+", feature_type=gm.ISOFORM_TYPE)
    objects = [
        gm.TranscriptRegion(RecordType.CDS, 1, 2, "chr1", "+"),
        Exon(1, 2, "chr1", "+"),
        isoform,
        gm.Transcript("TX1", 1, 20, "chr1", "+", feature_type=RecordType.MRNA),
        Gene("GENE1", None, 1, 20, "chr1", "+"),
        gm.PseudoGene("PSEUDO1", None, 1, 20, "chr1", "+"),
    ]
    for item in objects:
        assert not hasattr(item, "__dict__")


def test_pseudogene_reuses_gene_initialization_shape() -> None:
    pseudo = gm.PseudoGene("PSEUDO1", None, 1, 20, "chr1", "+")

    assert pseudo.feature_type is RecordType.PSEUDOGENE
    assert isinstance(pseudo.transcripts, dict)
    assert isinstance(pseudo.exons, list)
    assert isinstance(pseudo.cds, list)


def test_gene_default_attributes_are_not_shared() -> None:
    gene_one = Gene("GENE1", None, 1, 10, "chr1", "+")
    gene_two = Gene("GENE2", None, 20, 30, "chr1", "+")

    gene_one.attributes["tag"] = "one"

    assert "tag" not in gene_two.attributes


def test_isoform_constructor_does_not_contain_hardcoded_identifier_traps() -> None:
    isoform = gm.Transcript(
        "ENSG00000149256", 10, 20, "chr1", "+", feature_type=gm.ISOFORM_TYPE
    )

    assert isoform.id == "ENSG00000149256"


def test_basefeature_is_unhashable_because_coordinates_are_mutable() -> None:
    left = gm.BaseFeature("feature_a", 1, 10, "chr1", "+")
    with pytest.raises(TypeError, match="unhashable"):
        hash(left)


def test_load_gene_model_delegates_to_parser_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from SpliceGrapher.formats.parsers import gene_model_gff as parser_boundary

    model = GeneModel()
    records = ["chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1"]
    calls: dict[str, object] = {}

    def _fake_loader(
        model_arg: GeneModel, gff_records_arg: gm.GffRecordSource, **kwargs: object
    ) -> None:
        calls["model"] = model_arg
        calls["records"] = gff_records_arg
        calls["kwargs"] = kwargs

    monkeypatch.setattr(parser_boundary, "load_gene_model_records", _fake_loader)

    model.load_gene_model(records, verbose=True)

    assert calls["model"] is model
    assert calls["records"] == records
    assert calls["kwargs"] == {
        "require_notes": False,
        "chromosomes": None,
        "verbose": True,
        "ignore_errors": False,
    }


def test_load_gene_model_delegates_to_repository(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    model = GeneModel()
    records = ["chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1"]
    calls: dict[str, object] = {}

    def _fake_repo_load(
        model_arg: GeneModel,
        gff_records_arg: gm.GffRecordSource,
        *,
        require_notes: bool,
        chromosomes: Sequence[str] | str | None,
        verbose: bool,
        ignore_errors: bool,
    ) -> None:
        calls["model"] = model_arg
        calls["records"] = gff_records_arg
        calls["kwargs"] = {
            "require_notes": require_notes,
            "chromosomes": chromosomes,
            "verbose": verbose,
            "ignore_errors": ignore_errors,
        }

    monkeypatch.setattr(gm.GeneModelRepository, "load", _fake_repo_load)

    model.load_gene_model(records, verbose=True)

    assert calls["model"] is model
    assert calls["records"] == records
    assert calls["kwargs"] == {
        "require_notes": False,
        "chromosomes": None,
        "verbose": True,
        "ignore_errors": False,
    }


def test_load_gene_model_populates_mrna_parent_lookup() -> None:
    model = GeneModel()
    records = [
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=TX1;Parent=GENE1",
    ]

    model.load_gene_model(records)

    parent = model.get_mrna_parent("chr1", "tx1")
    assert parent is not None
    assert parent.id == "GENE1"


def test_write_gff_delegates_to_writer_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import io

    model = GeneModel()
    out_stream = io.StringIO()
    calls: dict[str, object] = {}

    def _fake_write_gff(
        model_arg: GeneModel,
        gff_path_arg: str | io.StringIO,
        *,
        gene_filter: gm.GeneFilter,
        gene_set: set[str] | list[str] | tuple[str, ...] | None,
        verbose: bool,
    ) -> None:
        calls["model"] = model_arg
        calls["path"] = gff_path_arg
        calls["gene_filter"] = gene_filter
        calls["gene_set"] = gene_set
        calls["verbose"] = verbose

    monkeypatch.setattr(gm, "write_gene_model_gff", _fake_write_gff)

    model.write_gff(out_stream, verbose=True)

    assert calls["model"] is model
    assert calls["path"] is out_stream
    assert calls["gene_filter"] is gm.default_gene_filter
    assert calls["gene_set"] is None
    assert calls["verbose"] is True


def test_write_gtf_delegates_to_writer_boundary(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    import io

    model = GeneModel()
    out_stream = io.StringIO()
    calls: dict[str, object] = {}

    def _fake_write_gtf(
        model_arg: GeneModel,
        gtf_path_arg: str | io.StringIO,
        *,
        gene_filter: gm.GeneFilter,
        verbose: bool,
    ) -> None:
        calls["model"] = model_arg
        calls["path"] = gtf_path_arg
        calls["gene_filter"] = gene_filter
        calls["verbose"] = verbose

    monkeypatch.setattr(gm, "write_gene_model_gtf", _fake_write_gtf)

    model.write_gtf(out_stream, verbose=True)

    assert calls["model"] is model
    assert calls["path"] is out_stream
    assert calls["gene_filter"] is gm.default_gene_filter
    assert calls["verbose"] is True


def test_gene_model_init_signature_uses_explicit_keywords() -> None:
    signature = inspect.signature(GeneModel.__init__)
    assert "args" not in signature.parameters
    assert "require_notes" in signature.parameters
    assert "chromosomes" in signature.parameters
    assert "verbose" in signature.parameters
    assert "ignore_errors" in signature.parameters


def test_gene_string_representation_is_not_stale_after_coordinate_updates() -> None:
    gene = Gene("GENE1", None, 10, 20, "chr1", "+")

    original = str(gene)
    gene.minpos = 1
    gene.maxpos = 30
    updated = str(gene)

    assert original != updated
    assert "1-30" in updated


def test_query_method_signatures_use_explicit_verbose_keywords() -> None:
    get_all_genes_sig = inspect.signature(GeneModel.get_all_genes)
    get_gene_records_sig = inspect.signature(GeneModel.get_gene_records)
    isoform_dict_sig = inspect.signature(GeneModel.isoform_dict)

    assert "args" not in get_all_genes_sig.parameters
    assert "verbose" in get_all_genes_sig.parameters
    assert "args" not in get_gene_records_sig.parameters
    assert "verbose" in get_gene_records_sig.parameters
    assert "args" not in isoform_dict_sig.parameters
    assert "verbose" in isoform_dict_sig.parameters


def test_get_parent_uses_strict_lookup_without_delimiter_guessing() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1")
    model.add_gene(gene)

    assert model.get_parent("GENE1", "chr1") is gene
    assert model.get_parent("GENE1.1", "chr1") is None


def test_load_gene_model_resolves_delimited_parent_ids_in_parser_layer() -> None:
    model = GeneModel()
    records = [
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsrc\texon\t10\t20\t.\t+\t.\tParent=GENE1.1;transcript_id=TX1",
    ]

    model.load_gene_model(records)

    gene = model.get_gene("chr1", "GENE1")
    assert gene is not None
    assert len(gene.exons) == 1
    assert "TX1" in gene.transcripts


def test_get_all_genes_keyword_is_supported() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    model.add_gene(Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1"))

    genes = model.get_all_genes(gene_filter=gm.default_gene_filter)

    assert [g.id for g in genes] == ["GENE1"]


def test_get_gene_by_name_is_case_insensitive_for_mixed_case_ids() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("At1g01160", None, 10, 20, "chr1", "+", name="At1g01160")
    model.add_gene(gene)

    assert model.get_gene_by_name("AT1G01160") is gene
    assert model.get_gene_by_name("at1g01160") is gene


def test_add_gene_rejects_case_only_duplicate_ids_within_chromosome() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    model.add_gene(Gene("At1g01160", None, 10, 20, "chr1", "+", name="At1g01160"))

    with pytest.raises(ValueError, match="already stored"):
        model.add_gene(Gene("AT1G01160", None, 30, 40, "chr1", "+", name="AT1G01160"))


def test_get_chromosome_normalizes_lookup_key_to_lowercase() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "Chr1")

    chromosome = model.get_chromosome("CHR1")

    assert chromosome is not None
    assert chromosome.name == "Chr1"


def test_add_gene_initializes_missing_chromosome_bucket_with_normalized_key() -> None:
    model = GeneModel()
    gene = Gene("GENE1", None, 10, 20, "Chr2", "+", name="GENE1")

    model.add_gene(gene)

    assert model.get_gene("chr2", "GENE1") is gene
    assert model.get_gene("CHR2", "GENE1") is gene


def test_get_parent_keywords_are_supported() -> None:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1")
    model.add_gene(gene)

    assert (
        model.get_parent("GENE1", "chr1", search_genes=True, search_mrna=False) is gene
    )
