from __future__ import annotations

from SpliceGrapher.core.enums import JunctionCode, ShortReadCode
from SpliceGrapher.shared import ShortRead as shortread


def test_shortread_codes_are_enum_backed() -> None:
    assert shortread.CHROM_CODE is ShortReadCode.CHROM
    assert shortread.DEPTH_CODE is ShortReadCode.DEPTH
    assert shortread.READ_CODE is ShortReadCode.READ
    assert shortread.JCT_CODE is ShortReadCode.JUNCTION
    assert shortread.SPLICE_CODE is ShortReadCode.SPLICE

    assert shortread.KNOWN_JCT is JunctionCode.KNOWN
    assert shortread.UNKNOWN_JCT is JunctionCode.UNKNOWN
    assert shortread.PREDICTED_JCT is JunctionCode.PREDICTED
    assert shortread.UNLABELED_JCT is JunctionCode.UNLABELED


def test_shortread_code_collections_are_enum_backed() -> None:
    assert all(isinstance(code, ShortReadCode) for code in shortread.KNOWN_CODES)
    assert all(isinstance(code, ShortReadCode) for code in shortread.DEPTH_CODES)
    assert all(isinstance(code, JunctionCode) for code in shortread.SJ_CODES)
