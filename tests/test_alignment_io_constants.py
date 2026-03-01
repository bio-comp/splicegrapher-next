from __future__ import annotations

from SpliceGrapher.core.enums import SamHeaderLine, SamHeaderTag
from SpliceGrapher.formats import alignment_io


def test_alignment_io_header_tags_are_enum_backed() -> None:
    assert alignment_io.HEADER_HD_TAG is SamHeaderTag.HD
    assert alignment_io.HEADER_LN_TAG is SamHeaderTag.LN
    assert alignment_io.HEADER_SN_TAG is SamHeaderTag.SN
    assert alignment_io.HEADER_SO_TAG is SamHeaderTag.SO
    assert alignment_io.HEADER_SQ_TAG is SamHeaderTag.SQ
    assert alignment_io.HEADER_VN_TAG is SamHeaderTag.VN
    assert alignment_io.HEADER_SQ_LINE is SamHeaderLine.SQ
