from __future__ import annotations

from SpliceGrapher.shared.format_utils import (
    comma_format,
    dict_string,
    list_string,
    substring_after,
    substring_before,
    substring_between,
    time_string,
    timestamp,
    to_numeric,
)


def test_comma_format_variants() -> None:
    assert comma_format("1234567") == "1,234,567"
    assert comma_format(8900) == "8,900"


def test_dict_and_list_string() -> None:
    mapping = {"a": 1, "b": 2}
    assert dict_string(mapping, delim="|") == "a -> 1|b -> 2"
    assert dict_string(mapping) == "a -> 1,b -> 2"

    values = ["x", 2, 3.5]
    assert list_string(values, delim=";") == "x;2;3.5"
    assert list_string(values) == "x,2,3.5"


def test_substring_helpers() -> None:
    text = "prefix[tag]middle[/tag]suffix"

    assert substring_after(text, "prefix") == "[tag]middle[/tag]suffix"
    assert substring_before(text, "suffix") == "prefix[tag]middle[/tag]"
    assert substring_between(text, "[tag]", "[/tag]") == "middle"

    assert substring_after(text, "missing") is None
    assert substring_before(text, "missing") is None
    assert substring_between(text, "[tag]", "missing") is None


def test_time_helpers() -> None:
    stamp = timestamp("%Y")
    assert len(stamp) == 4

    message = time_string("hello", format_string="%H:%M:%S")
    assert message.endswith(" hello")

    newline_message = time_string("world", format_string="%H:%M:%S", trailing_newline=True)
    assert newline_message.endswith(" world\n")


def test_to_numeric() -> None:
    assert to_numeric("10") == 10
    assert to_numeric("12.5") == 12.5
    assert to_numeric("abc") == "abc"
