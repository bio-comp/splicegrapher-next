from __future__ import annotations

import pytest

from SpliceGrapher.shared import collection_utils
from SpliceGrapher.shared.collection_utils import as_list, as_set, binary_search


def test_as_list_handles_supported_inputs() -> None:
    assert as_list("a,b,c") == ["a", "b", "c"]
    assert as_list(["a", "b"]) == ["a", "b"]
    assert as_list(("a", "b")) == ["a", "b"]
    assert sorted(as_list({"a", "b"})) == ["a", "b"]


def test_as_set_handles_supported_inputs() -> None:
    assert as_set("a,b,a") == {"a", "b"}
    assert as_set(["a", "b", "a"]) == {"a", "b"}
    assert as_set(("a", "b", "a")) == {"a", "b"}
    assert as_set({"a", "b"}) == {"a", "b"}


def test_as_list_and_as_set_reject_invalid_inputs() -> None:
    with pytest.raises(TypeError):
        as_list(123)  # type: ignore[arg-type]
    with pytest.raises(TypeError):
        as_set(123)  # type: ignore[arg-type]


def test_binary_search_returns_raw_insertion_index() -> None:
    sequence = [10, 20, 30]
    assert binary_search(sequence, 5) == 0
    assert binary_search(sequence, 20) == 1
    assert binary_search(sequence, 25) == 2
    assert binary_search(sequence, 40) == 3


def test_binary_search_supports_key() -> None:
    sequence = [("a", 10), ("b", 20), ("c", 30)]
    assert binary_search(sequence, 25, key=lambda item: item[1]) == 2


def test_binary_search_rejects_empty_sequence() -> None:
    with pytest.raises(ValueError, match="empty"):
        binary_search([], 1)


def test_legacy_wrapper_names_removed() -> None:
    assert not hasattr(collection_utils, "asList")
    assert not hasattr(collection_utils, "asSet")
    assert not hasattr(collection_utils, "bsearch")
