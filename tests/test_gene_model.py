from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pytest

from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.formats import GeneModel as gm
from SpliceGrapher.formats.GeneModel import Exon, Gene, GeneModel, feature_search, featureSearch


@dataclass
class _DummyGene:
    minpos: int
    maxpos: int
    strand: str = "+"

    def contains(self, loc: int, strand: str) -> bool:
        return strand == self.strand and self.minpos <= loc <= self.maxpos


def test_feature_search_supports_recursive_midpoint_indexing() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert featureSearch(features, query) is features[1]


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

    assert gm.featureOverlaps(left, right)
    assert gm.featureContains(left, nested)
    assert not gm.featureContains(left, right)


def test_feature_search_returns_preceding_feature_when_not_contained() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(32, 33, "chr1", "+")

    assert feature_search(features, query) is features[1]


def test_binary_search_genes_handles_large_gene_lists() -> None:
    model = GeneModel(None)
    genes = [_DummyGene(i * 10, (i * 10) + 9) for i in range(30)]

    low_gene, high_gene = model.binary_search_genes(genes, 155, 0, len(genes) - 1)

    assert low_gene is genes[15]
    assert high_gene is genes[15]


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


def test_legacy_camel_case_method_aliases_remain_available() -> None:
    model = GeneModel(None)
    assert model.addChromosome is not None
    assert model.loadGeneModel is not None
    assert model.makeSortedModel is not None
    assert model.writeGFF is not None


def test_gene_model_record_type_collections_are_enum_backed() -> None:
    assert all(isinstance(item, RecordType) for item in gm.KNOWN_RECTYPES)
    assert all(isinstance(item, RecordType) for item in gm.IGNORE_RECTYPES)
    assert all(isinstance(item, RecordType) for item in gm.CDS_TYPES)


def test_gene_model_alias_mapping_uses_recordtype_domain() -> None:
    assert gm.RECTYPE_MAP[RecordType.PREDICTED_GENE] is RecordType.GENE
    assert gm.RECTYPE_MAP[RecordType.CDS_PREDICTED] is RecordType.CDS


def test_gene_model_valid_strands_uses_enum_domain() -> None:
    assert gm.VALID_STRANDS == set(Strand)


def test_exon_default_attributes_are_not_shared() -> None:
    exon_one = Exon(1, 10, "chr1", "+")
    exon_two = Exon(20, 30, "chr1", "+")

    exon_one.attributes["tag"] = "one"

    assert "tag" not in exon_two.attributes


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
    model = GeneModel(None)
    records = ["chr1\tsrc\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1"]
    calls: dict[str, object] = {}

    from SpliceGrapher.formats.parsers import gene_model_gff

    def _fake_loader(
        model_arg: GeneModel, gff_records_arg: gm.GffRecordSource, **kwargs: object
    ) -> None:
        calls["model"] = model_arg
        calls["records"] = gff_records_arg
        calls["kwargs"] = kwargs

    monkeypatch.setattr(gene_model_gff, "load_gene_model_records", _fake_loader)

    model.load_gene_model(records, verbose=True)

    assert calls["model"] is model
    assert calls["records"] == records
    assert calls["kwargs"] == {"verbose": True}


def test_write_gff_delegates_to_writer_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    import io

    model = GeneModel(None)
    out_stream = io.StringIO()
    calls: dict[str, object] = {}

    from SpliceGrapher.formats.writers import gene_model

    def _fake_write_gff(
        model_arg: GeneModel, gff_path_arg: str | io.StringIO, **kwargs: object
    ) -> None:
        calls["model"] = model_arg
        calls["path"] = gff_path_arg
        calls["kwargs"] = kwargs

    monkeypatch.setattr(gene_model, "write_gff", _fake_write_gff)

    model.write_gff(out_stream, verbose=True)

    assert calls["model"] is model
    assert calls["path"] is out_stream
    assert calls["kwargs"] == {"verbose": True}


def test_write_gtf_delegates_to_writer_boundary(monkeypatch: pytest.MonkeyPatch) -> None:
    import io

    model = GeneModel(None)
    out_stream = io.StringIO()
    calls: dict[str, object] = {}

    from SpliceGrapher.formats.writers import gene_model

    def _fake_write_gtf(
        model_arg: GeneModel,
        gtf_path_arg: str | io.StringIO,
        geneFilter: gm.GeneFilter = gm.defaultGeneFilter,
        **kwargs: object,
    ) -> None:
        calls["model"] = model_arg
        calls["path"] = gtf_path_arg
        calls["gene_filter"] = geneFilter
        calls["kwargs"] = kwargs

    monkeypatch.setattr(gene_model, "write_gtf", _fake_write_gtf)

    model.write_gtf(out_stream, verbose=True)

    assert calls["model"] is model
    assert calls["path"] is out_stream
    assert calls["gene_filter"] is gm.defaultGeneFilter
    assert calls["kwargs"] == {"verbose": True}
