from __future__ import annotations

import importlib
import inspect
from pathlib import Path

import pytest

from SpliceGrapher.formats.alignment_io import collect_alignment_data
from tests.helpers.alignment_fixture_builder import build_alignment_fixture

alignment_api = importlib.import_module("SpliceGrapher.formats.alignment_io.api")
alignment_collect = importlib.import_module("SpliceGrapher.formats.alignment_io.collect")
alignment_depths = importlib.import_module("SpliceGrapher.formats.alignment_io.depths")
alignment_sources = importlib.import_module("SpliceGrapher.formats.alignment_io.sources")


def _alignment_package_source() -> str:
    package_dir = Path(__file__).resolve().parents[1] / "SpliceGrapher" / "formats" / "alignment_io"
    return "\n".join(path.read_text(encoding="utf-8") for path in sorted(package_dir.glob("*.py")))


def test_alignment_io_no_longer_uses_process_utils_getattribute() -> None:
    source = _alignment_package_source()
    assert "from SpliceGrapher.shared.process_utils import getAttribute" not in source
    assert "getAttribute(" not in source


def test_alignment_io_removes_manual_record_parser_and_list_growth_antipattern() -> None:
    source = _alignment_package_source()
    assert "class AlignmentRecord" not in source
    assert "+= [0] * len(" not in source
    assert ".tolist()" not in source


@pytest.mark.parametrize(
    ("module", "func_name"),
    [
        (alignment_sources, "_open_alignment_file"),
        (alignment_collect, "_collect_pysam_data"),
        (alignment_depths, "calculate_gene_depths"),
        (alignment_api, "read_alignment_spans"),
        (alignment_api, "read_alignment_depths"),
        (alignment_api, "read_alignment_headers"),
        (alignment_api, "read_alignment_chromosome_info"),
        (alignment_api, "read_alignment_junctions"),
        (alignment_api, "collect_alignment_data"),
    ],
)
def test_alignment_io_hot_path_signatures_are_explicit(module: object, func_name: str) -> None:
    signature = inspect.signature(getattr(module, func_name))
    assert all(
        param.kind != inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values()
    )


def test_get_sam_read_data_rejects_unknown_keyword_argument(tmp_path: Path) -> None:
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)
    with pytest.raises(TypeError):
        collect_alignment_data(  # type: ignore[call-overload]
            str(fixture.bam), nonsense_option=True
        )


def test_get_sam_read_data_depths_fallback_preserves_three_tuple_when_alignments_requested(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.depths"
    sample_path.write_text("placeholder\n", encoding="utf-8")

    fake_depths = {"chr1": [0, 1, 2]}
    fake_junctions: dict[str, list[object]] = {"chr1": []}

    monkeypatch.setattr(alignment_api, "_is_depths_source", lambda source: True)
    monkeypatch.setattr(
        alignment_api,
        "read_depths",
        lambda *args, **kwargs: (fake_depths, fake_junctions),
    )

    depths, junctions, alignments = collect_alignment_data(
        str(sample_path), include_alignments=True
    )

    assert list(depths["chr1"]) == [0, 1, 2]
    assert junctions == fake_junctions
    assert alignments == {}
