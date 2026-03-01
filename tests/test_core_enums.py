from __future__ import annotations

import pytest

from SpliceGrapher.core import enum_coercion
from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType, Strand


def test_strand_values_are_stable() -> None:
    assert Strand.PLUS.value == "+"
    assert Strand.MINUS.value == "-"
    assert Strand.UNKNOWN.value == "."


def test_coerce_enum_accepts_symbol_for_strand() -> None:
    assert coerce_enum("+", Strand, field="strand") is Strand.PLUS


def test_coerce_enum_rejects_invalid_strand() -> None:
    with pytest.raises(ValueError):
        coerce_enum("?", Strand, field="strand")


def test_coerce_enum_rejects_unknown_record_type() -> None:
    with pytest.raises(ValueError):
        coerce_enum("totally_unknown_type", RecordType, field="record_type")


def test_wrapper_helpers_removed() -> None:
    assert not hasattr(enum_coercion, "coerce_strand")
    assert not hasattr(enum_coercion, "coerce_record_type")
    assert not hasattr(enum_coercion, "coerce_node_disposition")
    assert not hasattr(enum_coercion, "coerce_edge_type")
    assert not hasattr(enum_coercion, "coerce_attr_key")
