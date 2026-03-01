from __future__ import annotations

import pytest

from SpliceGrapher.core.enum_coercion import coerce_record_type, coerce_strand
from SpliceGrapher.core.enums import Strand


def test_strand_values_are_stable() -> None:
    assert Strand.PLUS.value == "+"
    assert Strand.MINUS.value == "-"
    assert Strand.UNKNOWN.value == "."


def test_coerce_strand_accepts_symbol() -> None:
    assert coerce_strand("+") is Strand.PLUS


def test_coerce_strand_rejects_invalid() -> None:
    with pytest.raises(ValueError):
        coerce_strand("?")


def test_coerce_record_type_rejects_unknown() -> None:
    with pytest.raises(ValueError):
        coerce_record_type("totally_unknown_type")
