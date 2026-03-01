from __future__ import annotations

from SpliceGrapher.core.enums import Strand


def test_strand_values_are_stable() -> None:
    assert Strand.PLUS.value == "+"
    assert Strand.MINUS.value == "-"
    assert Strand.UNKNOWN.value == "."
