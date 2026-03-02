from __future__ import annotations

from SpliceGrapher.core.enums import SamHeaderTag
from SpliceGrapher.formats import alignment_io


def test_alignment_io_uses_direct_header_enums_without_proxy_exports() -> None:
    assert not hasattr(alignment_io, "HEADER_HD_TAG")
    assert not hasattr(alignment_io, "HEADER_LN_TAG")
    assert not hasattr(alignment_io, "HEADER_SN_TAG")
    assert not hasattr(alignment_io, "HEADER_SO_TAG")
    assert not hasattr(alignment_io, "HEADER_SQ_TAG")
    assert not hasattr(alignment_io, "HEADER_VN_TAG")
    assert not hasattr(alignment_io, "HEADER_SQ_LINE")

    class _HeaderStub:
        header = {
            SamHeaderTag.HD: {
                SamHeaderTag.VN: "1.0",
                SamHeaderTag.SO: "coordinate",
            },
            SamHeaderTag.SQ: [
                {SamHeaderTag.SN: "chr1", SamHeaderTag.LN: 1000},
            ],
        }

    rendered = alignment_io.pysamHeaders(_HeaderStub())
    assert rendered[0].startswith("@HD")
    assert rendered[1].startswith("@SQ")


def test_alignment_io_does_not_export_maxint_alias() -> None:
    assert not hasattr(alignment_io, "MAXINT")
