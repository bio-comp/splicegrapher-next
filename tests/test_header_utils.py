from __future__ import annotations

import pytest

from SpliceGrapher.shared.header_utils import (
    process_fasta_header,
    process_fastq_header,
    process_labeled_fasta_header,
)


def test_process_fastq_header_returns_first_token() -> None:
    assert process_fastq_header("chr1 read=abc") == "chr1"


def test_process_fasta_header_returns_first_token() -> None:
    assert process_fasta_header("chr2 transcript=xyz") == "chr2"


def test_process_fasta_like_header_rejects_empty_values() -> None:
    with pytest.raises(ValueError, match="Header is empty"):
        process_fastq_header("   ")

    with pytest.raises(ValueError, match="Header is empty"):
        process_fasta_header("")


def test_process_labeled_fasta_header_parses_label_with_extra_tokens() -> None:
    seq_id, label = process_labeled_fasta_header("seqA label=retained source=training")
    assert seq_id == "seqA"
    assert label == "retained"


def test_process_labeled_fasta_header_requires_label_token() -> None:
    with pytest.raises(ValueError, match="label="):
        process_labeled_fasta_header("seqA source=training")
