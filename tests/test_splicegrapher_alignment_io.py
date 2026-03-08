from pathlib import Path

import numpy

from SpliceGrapher.formats.alignment_io import (
    collect_alignment_data,
    read_alignment_depths,
    read_alignment_headers,
    read_alignment_junctions,
    read_alignment_spans,
)
from tests.helpers.alignment_fixture_builder import build_alignment_fixture


def _junction_signature(junctions_by_chrom):
    """Build a stable signature for comparing junction collections."""
    signature = {}
    for chrom, junctions in junctions_by_chrom.items():
        signature[chrom] = sorted((j.donor(), j.acceptor(), j.strand, j.count) for j in junctions)
    return signature


def _depth_maps_match(
    left_depths: dict[str, numpy.ndarray],
    right_depths: dict[str, numpy.ndarray],
) -> bool:
    if left_depths.keys() != right_depths.keys():
        return False
    return all(numpy.array_equal(left_depths[chrom], right_depths[chrom]) for chrom in left_depths)


def test_get_sam_read_data_matches_between_sam_and_bam(tmp_path: Path):
    """SAM and BAM inputs should produce equivalent depth/junction summaries."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=24)

    sam_depths, sam_junctions = collect_alignment_data(str(fixture.sam))
    bam_depths, bam_junctions = collect_alignment_data(str(fixture.bam))

    assert _depth_maps_match(sam_depths, bam_depths)
    assert _junction_signature(sam_junctions) == _junction_signature(bam_junctions)


def test_get_sam_read_data_supports_cram_with_reference(tmp_path: Path):
    """CRAM reads should load when an explicit reference FASTA is provided."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=24)

    bam_depths, bam_junctions = collect_alignment_data(str(fixture.bam))
    cram_depths, cram_junctions = collect_alignment_data(
        str(fixture.cram),
        reference_fasta=str(fixture.reference_fasta),
    )

    assert _depth_maps_match(cram_depths, bam_depths)
    assert _junction_signature(cram_junctions) == _junction_signature(bam_junctions)


def test_get_sam_read_data_cram_with_bad_reference_is_actionable(tmp_path: Path):
    """CRAM decoding should either succeed or fail with a clear reference-related error."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)
    bad_reference = tmp_path / "missing.fa"

    try:
        depths, junctions = collect_alignment_data(
            str(fixture.cram), reference_fasta=str(bad_reference)
        )
        assert depths
        assert fixture.chrom in depths
        assert fixture.chrom in junctions
    except Exception as exc:
        message = str(exc).lower()
        assert "reference" in message or "cram" in message or "no such file" in message


def test_get_sam_depths_and_junctions_respect_chromosome_filter(tmp_path: Path):
    """Depth and junction helpers should respect the chromosome filter argument."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=16)

    depths = read_alignment_depths(str(fixture.bam), chromosomes=[fixture.chrom])
    junctions = read_alignment_junctions(str(fixture.bam), chromosomes=[fixture.chrom], minjct=1)

    assert set(depths.keys()) == {fixture.chrom}
    assert set(junctions.keys()) == {fixture.chrom}
    assert junctions[fixture.chrom]


def test_get_sam_read_data_alignments_true_returns_alignment_map(tmp_path: Path) -> None:
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)

    depths, junctions, alignments = collect_alignment_data(
        str(fixture.bam), include_alignments=True
    )

    assert fixture.chrom in depths
    assert fixture.chrom in junctions
    assert fixture.chrom in alignments
    assert alignments[fixture.chrom]


def test_get_sam_alignments_matches_read_data_alignment_map(tmp_path: Path) -> None:
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)

    _, _, alignments = collect_alignment_data(str(fixture.bam), include_alignments=True)
    direct_alignments = read_alignment_spans(str(fixture.bam))

    assert direct_alignments == alignments


def test_get_sam_headers_reads_bam_and_sam_headers(tmp_path: Path):
    """Header helper should work for both SAM text and BAM binary inputs."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)

    sam_headers = read_alignment_headers(str(fixture.sam))
    bam_headers = read_alignment_headers(str(fixture.bam))

    assert any(line.startswith("@HD") for line in sam_headers)
    assert any(line.startswith("@SQ") for line in sam_headers)
    assert any(line.startswith("@HD") for line in bam_headers)
    assert any(line.startswith("@SQ") for line in bam_headers)
