from __future__ import annotations

import importlib

from SpliceGrapher.formats import alignment_io

alignment_api = importlib.import_module("SpliceGrapher.formats.alignment_io.api")
alignment_collect = importlib.import_module("SpliceGrapher.formats.alignment_io.collect")
alignment_sources = importlib.import_module("SpliceGrapher.formats.alignment_io.sources")
alignment_types = importlib.import_module("SpliceGrapher.formats.alignment_io.types")


def test_alignment_io_package_reexports_public_api_and_internal_layout() -> None:
    assert alignment_io.collect_alignment_data is alignment_api.collect_alignment_data
    assert alignment_io.read_alignment_depths is alignment_api.read_alignment_depths
    assert hasattr(alignment_collect, "_collect_pysam_data")
    assert hasattr(alignment_sources, "_open_alignment_file")
    assert hasattr(alignment_types, "AlignmentSource")
    assert not hasattr(alignment_collect, "_is_depths_source")
    assert not hasattr(alignment_collect, "_collect_depths_source_data")
    assert hasattr(alignment_api, "_is_depths_source")
    assert hasattr(alignment_api, "_collect_depths_source_data")


def test_alignment_io_uses_direct_header_enums_without_proxy_exports() -> None:
    assert not hasattr(alignment_io, "HEADER_HD_TAG")
    assert not hasattr(alignment_io, "HEADER_LN_TAG")
    assert not hasattr(alignment_io, "HEADER_SN_TAG")
    assert not hasattr(alignment_io, "HEADER_SO_TAG")
    assert not hasattr(alignment_io, "HEADER_SQ_TAG")
    assert not hasattr(alignment_io, "HEADER_VN_TAG")
    assert not hasattr(alignment_io, "HEADER_SQ_LINE")


def test_alignment_io_does_not_export_maxint_alias() -> None:
    assert not hasattr(alignment_io, "MAXINT")
