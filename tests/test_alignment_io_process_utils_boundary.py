from __future__ import annotations

import inspect
from pathlib import Path

import pytest

from SpliceGrapher.formats import alignment_io
from SpliceGrapher.formats.alignment_io import getSamReadData
from tests.helpers.alignment_fixture_builder import build_alignment_fixture


def test_alignment_io_no_longer_uses_process_utils_getattribute() -> None:
    source = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "alignment_io.py"
    ).read_text(encoding="utf-8")
    assert "from SpliceGrapher.shared.process_utils import getAttribute" not in source
    assert "getAttribute(" not in source


def test_alignment_io_removes_manual_record_parser_and_list_growth_antipattern() -> None:
    source = (
        Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "alignment_io.py"
    ).read_text(encoding="utf-8")
    assert "class AlignmentRecord" not in source
    assert "+= [0] * len(" not in source


@pytest.mark.parametrize(
    "func_name",
    [
        "_open_alignment_file",
        "_collect_pysam_data",
        "pysamReadDepths",
        "getSamAlignments",
        "getSamDepths",
        "getSamHeaders",
        "getSamHeaderInfo",
        "getSamJunctions",
        "getSamReadData",
    ],
)
def test_alignment_io_hot_path_signatures_are_explicit(func_name: str) -> None:
    signature = inspect.signature(getattr(alignment_io, func_name))
    assert all(
        param.kind != inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()
    )


def test_get_sam_read_data_rejects_unknown_keyword_argument(tmp_path: Path) -> None:
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)
    with pytest.raises(TypeError):
        getSamReadData(str(fixture.bam), nonsense_option=True)
