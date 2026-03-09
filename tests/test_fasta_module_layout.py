from __future__ import annotations

from SpliceGrapher.formats import fasta


def test_fasta_is_package_with_modern_exports() -> None:
    assert hasattr(fasta, "__path__")
    assert hasattr(fasta, "FastaIterator")
    assert hasattr(fasta, "FastaSlice")
    assert hasattr(fasta, "truncate_sequences")
    assert not hasattr(fasta, "fasta_itr")
    assert not hasattr(fasta, "fasta_slice")
    assert not hasattr(fasta, "truncateSequences")
