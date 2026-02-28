from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from SpliceGrapher.formats import fasta


def _write_fasta(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    return path


def test_fasta_iterator_accepts_pathlib_paths(tmp_path: Path) -> None:
    fasta_path = _write_fasta(
        tmp_path / "example.fa",
        ">chr1\nACGT\nAC\n>chr2\nTTAA\n",
    )

    records = list(fasta.fasta_itr(fasta_path))

    assert [(record.header, record.sequence) for record in records] == [
        ("chr1", "ACGTAC"),
        ("chr2", "TTAA"),
    ]
    assert fasta.fasta_count(fasta_path) == 2


def test_fasta_iterator_accepts_open_file_handles(tmp_path: Path) -> None:
    fasta_path = _write_fasta(tmp_path / "handle.fa", ">seqA\nAACCGG\n")

    with fasta_path.open("r", encoding="utf-8") as handle:
        records = list(fasta.fasta_itr(handle))

    assert len(records) == 1
    assert records[0].header == "seqA"
    assert records[0].sequence == "AACCGG"


def test_fasta_iterator_raises_on_blank_line(tmp_path: Path) -> None:
    fasta_path = _write_fasta(tmp_path / "bad.fa", ">seqA\nAA\n\nTT\n")

    with pytest.raises(fasta.MalformedInput):
        list(fasta.fasta_itr(fasta_path))


def test_get_sequence_preserves_long_headers(tmp_path: Path) -> None:
    fasta_path = _write_fasta(
        tmp_path / "headers.fa",
        ">chr1 transcript A\nAAAA\n>chr2 transcript B\nCCCC\n",
    )

    record = fasta.get_sequence(fasta_path, "chr1 transcript A")
    by_prefix = fasta.fasta_get_by_name(
        iter(fasta.fasta_itr(fasta_path)),
        "chr1",
        byLength=True,
    )

    assert record is not None
    assert record.header == "chr1 transcript A"
    assert record.sequence == "AAAA"
    assert by_prefix is not None
    assert by_prefix.header == "chr1 transcript A"


def test_fasta_iterator_reads_gzipped_input(tmp_path: Path) -> None:
    fasta_path = tmp_path / "compressed.fa.gz"
    with gzip.open(fasta_path, "wt", encoding="utf-8") as handle:
        handle.write(">gz_rec\nACGT\n")

    records = list(fasta.fasta_itr(str(fasta_path)))

    assert len(records) == 1
    assert records[0].header == "gz_rec"
    assert records[0].sequence == "ACGT"
