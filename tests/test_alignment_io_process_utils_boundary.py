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
    assert ".tolist()" not in source


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
        getattr(alignment_io, "getSamReadData")(str(fixture.bam), nonsense_option=True)


def test_get_sam_read_data_depths_fallback_preserves_three_tuple_when_alignments_requested(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.depths"
    sample_path.write_text("placeholder\n", encoding="utf-8")

    fake_depths = {"chr1": [0, 1, 2]}
    fake_junctions: dict[str, list[object]] = {"chr1": []}

    monkeypatch.setattr(alignment_io, "_is_depths_source", lambda source: True)
    monkeypatch.setattr(
        alignment_io,
        "read_depths",
        lambda *args, **kwargs: (fake_depths, fake_junctions),
    )

    depths, junctions, alignments = getSamReadData(str(sample_path), alignments=True)

    assert list(depths["chr1"]) == [0, 1, 2]
    assert junctions == fake_junctions
    assert alignments == {}
