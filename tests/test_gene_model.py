from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.formats import gene_model as gm
from SpliceGrapher.formats.gene_model import Exon, Gene, GeneModel, feature_search


def test_feature_search_supports_recursive_midpoint_indexing() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert feature_search(features, query) is features[1]


def test_mrna_sorted_exons_uses_explicit_minintron_keyword() -> None:
    transcript = gm.Mrna("tx1", 1, 20, "chr1", "+")
    transcript.add_cds(gm.CDS(1, 5, "chr1", "+"))
    transcript.add_cds(gm.CDS(7, 10, "chr1", "+"))

    merged = transcript.sorted_exons(minintron=3)
    unmerged = transcript.sorted_exons(minintron=2)

    assert len(merged) == 1
    assert len(unmerged) == 2

    with pytest.raises(TypeError):
        transcript.sorted_exons(min_intron=3)  # type: ignore[call-arg]


def test_feature_search_snake_case_api_returns_expected_match() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert feature_search(features, query) is features[1]


def test_feature_overlap_and_contains_match_interval_helper_semantics() -> None:
    left = Exon(10, 20, "chr1", "+")
    right = Exon(20, 30, "chr1", "+")
    nested = Exon(12, 18, "chr1", "+")

    assert gm.feature_overlaps(left, right)
    assert gm.feature_contains(left, nested)
    assert not gm.feature_contains(left, right)


def test_feature_search_returns_preceding_feature_when_not_contained() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(32, 33, "chr1", "+")

    assert feature_search(features, query) is features[1]


def test_get_gene_from_locations_uses_interval_index() -> None:
    model = GeneModel(None)
    model.add_chromosome(1, 1000, "chr1")
    gene = Gene("GENE_A", None, 150, 190, "chr1", "+", name="GENE_A")
    model.add_gene(gene)
    model.make_sorted_model()

    assert model.get_gene_from_locations("chr1", 155, 160, "+") is gene


def test_write_gff_writes_chromosome_records(tmp_path: Path) -> None:
    model = GeneModel(None)
    model.add_chromosome(1, 100, "chr1")
    model.add_gene(
        Gene(
            id="GENE1",
            note=None,
            start=10,
            end=20,
            chromosome="chr1",
            strand="+",
            name="GENE1",
            attr={},
        )
    )

    gff_path = tmp_path / "model.gff3"
    model.write_gff(str(gff_path))

    lines = gff_path.read_text(encoding="utf-8").splitlines()
    expected_prefix = "Chr1\tSpliceGrapher\tchromosome\t1\t100\t.\t.\t.\tID=Chr1;Name=Chr1"
    assert lines[0].startswith(expected_prefix)


def test_get_annotation_dict_ignores_malformed_tokens_and_splits_once() -> None:
    model = GeneModel(None)
    parsed = model.get_annotation_dict("ID=GENE1;badtoken;Name=A=B;Parent=P1;=junk;")

    assert parsed == {"ID": "GENE1", "Name": "A=B", "Parent": "P1"}


def test_load_gene_model_rejects_unknown_record_type() -> None:
    model = GeneModel(None)
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
    isoform = gm.Isoform("TX1", 1, 200, "chr1", "-")
    isoform.set_parent(gene)
    isoform.add_exon(Exon(80, 90, "chr1", "-"))
    isoform.add_exon(Exon(20, 30, "chr1", "-"))

    starts = [int(line.split("\t")[3]) for line in isoform.gtf_strings()]
    assert starts == [20, 80]


def test_mrna_gtf_strings_remain_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    transcript = gm.Mrna("TX1", 1, 200, "chr1", "-")
    transcript.parent = gene
    transcript.add_cds(gm.CDS(80, 90, "chr1", "-"))
    transcript.add_cds(gm.CDS(20, 30, "chr1", "-"))

    cds_starts = [int(line.split("\t")[3]) for line in transcript.gtf_strings()]
    assert cds_starts == [20, 80]


def test_gene_gtf_strings_common_path_remains_genomic_ascending_on_minus_strand() -> None:
    gene = Gene("GENE1", None, 1, 200, "chr1", "-", name="GENE1")
    isoform = gm.Isoform("TX1", 1, 200, "chr1", "-")
    transcript = gm.Mrna("TX1", 1, 200, "chr1", "-")
    gene.add_isoform(isoform)
    gene.add_mrna(transcript)

    isoform.add_exon(Exon(70, 75, "chr1", "-"))
    isoform.add_exon(Exon(20, 30, "chr1", "-"))
    transcript.add_cds(gm.CDS(80, 90, "chr1", "-"))
    transcript.add_cds(gm.CDS(25, 35, "chr1", "-"))

    feature_starts = [
        int(line.split("\t")[3])
        for line in gene.gtf_strings().splitlines()
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
    isoform = gm.Isoform("TX1", 1, 100, "chr1", "+")
    exon = Exon(10, 20, "chr1", "+")

    gene.add_isoform(isoform)
    isoform.add_exon(exon)

    assert not hasattr(exon, "parents")


def test_gene_default_attributes_are_not_shared() -> None:
    gene_one = Gene("GENE1", None, 1, 10, "chr1", "+")
    gene_two = Gene("GENE2", None, 20, 30, "chr1", "+")

    gene_one.attributes["tag"] = "one"

    assert "tag" not in gene_two.attributes


def test_isoform_constructor_does_not_contain_hardcoded_identifier_traps() -> None:
    isoform = gm.Isoform("ENSG00000149256", 10, 20, "chr1", "+")

    assert isoform.id == "ENSG00000149256"


def test_basefeature_hash_matches_equality_contract() -> None:
    left = gm.BaseFeature("feature_a", 1, 10, "chr1", "+")
    right = gm.BaseFeature("feature_b", 1, 10, "chr1", "+")

    assert left == right
    assert hash(left) == hash(right)


def test_load_gene_model_delegates_to_parser_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    from SpliceGrapher.formats.parsers import gene_model_gff as parser_boundary

    model = GeneModel(None)
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


def test_load_gene_model_populates_mrna_parent_lookup() -> None:
    model = GeneModel(None)
    records = [
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsrc\tmRNA\t1\t100\t.\t+\t.\tID=TX1;Parent=GENE1",
    ]

    model.load_gene_model(records)

    parent = model.get_mrna_parent("chr1", "tx1")
    assert parent is not None
    assert parent.id == "GENE1"


def test_write_gff_delegates_to_writer_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    import io

    model = GeneModel(None)
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


def test_write_gtf_delegates_to_writer_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    import io

    model = GeneModel(None)
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
    model = GeneModel(None)
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1")
    model.add_gene(gene)

    assert model.get_parent("GENE1", "chr1") is gene
    assert model.get_parent("GENE1.1", "chr1") is None


def test_load_gene_model_resolves_delimited_parent_ids_in_parser_layer() -> None:
    model = GeneModel(None)
    records = [
        "chr1\tsrc\tgene\t1\t100\t.\t+\t.\tID=GENE1;Name=GENE1",
        "chr1\tsrc\texon\t10\t20\t.\t+\t.\tParent=GENE1.1;transcript_id=TX1",
    ]

    model.load_gene_model(records)

    gene = model.get_gene("chr1", "GENE1")
    assert gene is not None
    assert len(gene.exons) == 1
    assert "TX1" in gene.isoforms


def test_get_all_genes_keyword_is_supported() -> None:
    model = GeneModel(None)
    model.add_chromosome(1, 100, "chr1")
    model.add_gene(Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1"))

    genes = model.get_all_genes(gene_filter=gm.default_gene_filter)

    assert [g.id for g in genes] == ["GENE1"]


def test_get_parent_keywords_are_supported() -> None:
    model = GeneModel(None)
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("GENE1", None, 10, 20, "chr1", "+", name="GENE1")
    model.add_gene(gene)

    assert model.get_parent("GENE1", "chr1", search_genes=True, search_mrna=False) is gene
