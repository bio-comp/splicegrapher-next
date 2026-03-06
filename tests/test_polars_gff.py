from __future__ import annotations

from pathlib import Path

import pytest

import SpliceGrapher.formats.polars_gff as polars_gff
from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.formats.gene_model import CDS, Exon, Gene, GeneModel, Transcript


def test_parse_gff_attributes_ignores_malformed_tokens_and_splits_once() -> None:
    parsed = polars_gff.parse_gff_attributes("ID=GENE1;badtoken;Name=A=B;Parent=P1;=junk;")
    assert parsed == {"ID": "GENE1", "Name": "A=B", "Parent": "P1"}


def test_iter_gff_records_skips_comments_and_malformed_rows() -> None:
    records = list(
        polars_gff.iter_gff_records(
            [
                "# comment\n",
                "chr1\tsource\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1\n",
                "malformed\trow\n",
                "chr1\tsource\texon\tbad\t20\t.\t+\t.\tID=EX1;Parent=TX1\n",
                "chr1\tsource\texon\t11\t20\t.\t+\t.\tID=EX2;Parent=TX1;Name=A=B\n",
            ],
            ignore_malformed=True,
        )
    )

    assert len(records) == 2
    assert records[0]["feature_id"] == "GENE1"
    assert records[1]["feature_id"] == "EX2"
    assert records[1]["name"] == "A=B"


def test_load_gff_rows_reads_from_path(tmp_path: Path) -> None:
    gff_path = tmp_path / "sample.gff3"
    gff_path.write_text(
        "# test\n"
        "chr1\tsource\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1\n"
        "chr1\tsource\texon\t11\t20\t.\t+\t.\tID=EX2;Parent=TX1\n",
        encoding="utf-8",
    )

    rows = polars_gff.load_gff_rows(gff_path)
    assert len(rows) == 2
    assert rows[0]["feature_id"] == "GENE1"
    assert rows[1]["parent_id"] == "TX1"


def test_load_gff_to_polars_raises_without_polars(monkeypatch, tmp_path: Path) -> None:
    gff_path = tmp_path / "sample.gff3"
    gff_path.write_text(
        "chr1\tsource\tgene\t1\t10\t.\t+\t.\tID=GENE1;Name=GENE1\n",
        encoding="utf-8",
    )

    def _raise_import_error(name: str):
        raise ModuleNotFoundError(name)

    monkeypatch.setattr(polars_gff.importlib, "import_module", _raise_import_error)

    with pytest.raises(polars_gff.PolarsNotInstalledError):
        polars_gff.load_gff_to_polars(gff_path)


def test_iter_gff_records_normalizes_record_type_and_strand_domains() -> None:
    records = list(
        polars_gff.iter_gff_records(
            [
                "chr1\tsource\tmRNA\t1\t10\t.\t+\t.\tID=TX1;Parent=GENE1\n",
            ]
        )
    )

    assert records[0]["type"] == RecordType.MRNA
    assert records[0]["strand"] == Strand.PLUS


def test_iter_gff_records_rejects_unknown_record_type() -> None:
    with pytest.raises(ValueError, match="line 1: .*record_type"):
        list(
            polars_gff.iter_gff_records(
                ["chr1\tsource\tunknown_type\t1\t10\t.\t+\t.\tID=TX1;Parent=GENE1\n"]
            )
        )


def test_iter_gff_records_rejects_unknown_strand() -> None:
    with pytest.raises(ValueError, match="line 1: .*strand"):
        list(
            polars_gff.iter_gff_records(
                ["chr1\tsource\tgene\t1\t10\t.\t?\t.\tID=GENE1;Name=GENE1\n"]
            )
        )


def _build_gene_model_for_extractor() -> GeneModel:
    model = GeneModel()
    model.add_chromosome(1, 100, "chr1")
    gene = Gene("GENE1", None, 10, 40, "chr1", "+", name="GENE1")
    transcript = Transcript("TX1", 10, 40, "chr1", "+")
    transcript.add_exon(Exon(10, 18, "chr1", "+"))
    transcript.add_exon(Exon(20, 40, "chr1", "+"))
    transcript.add_cds(CDS(12, 15, "chr1", "+"))
    transcript.add_cds(CDS(20, 30, "chr1", "+"))
    gene.add_transcript(transcript)
    model.add_gene(gene)
    return model


def test_iter_flattened_features_emits_gene_transcript_and_child_rows() -> None:
    rows = list(polars_gff._iter_flattened_features(_build_gene_model_for_extractor()))

    feature_types = {str(row["feature_type"]) for row in rows}
    assert {"GENE", "MRNA", "EXON", "CDS"}.issubset(feature_types)
    assert all(row["gene_id"] == "GENE1" for row in rows)


def test_extract_to_dataframe_raises_without_polars(monkeypatch) -> None:
    def _raise_import_error(name: str):
        raise ModuleNotFoundError(name)

    monkeypatch.setattr(polars_gff.importlib, "import_module", _raise_import_error)

    with pytest.raises(polars_gff.PolarsNotInstalledError):
        polars_gff.extract_to_dataframe(_build_gene_model_for_extractor())


def test_extract_to_dataframe_builds_expected_schema_when_polars_available() -> None:
    polars = pytest.importorskip("polars")
    dataframe = polars_gff.extract_to_dataframe(_build_gene_model_for_extractor())

    assert dataframe.height > 0
    assert dataframe.schema["chromosome"] == polars.Categorical
    assert dataframe.schema["start_pos"] == polars.Int64
