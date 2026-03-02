from __future__ import annotations

from pathlib import Path

import pytest

import SpliceGrapher.formats.polars_gff as polars_gff
from SpliceGrapher.core.enums import RecordType, Strand


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
