from pathlib import Path

from SpliceGrapher.formats.alignment_io import (
    getSamDepths,
    getSamJunctions,
    getSamReadData,
)
from tests.helpers.alignment_fixture_builder import build_alignment_fixture


def _junction_signature(junctions_by_chrom):
    """Build a stable signature for comparing junction collections."""
    signature = {}
    for chrom, junctions in junctions_by_chrom.items():
        signature[chrom] = sorted((j.donor(), j.acceptor(), j.strand, j.count) for j in junctions)
    return signature


def test_get_sam_read_data_matches_between_sam_and_bam(tmp_path: Path):
    """SAM and BAM inputs should produce equivalent depth/junction summaries."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=24)

    sam_depths, sam_junctions = getSamReadData(str(fixture.sam))
    bam_depths, bam_junctions = getSamReadData(str(fixture.bam))

    assert sam_depths.keys() == bam_depths.keys()
    chrom = fixture.chrom
    assert sam_depths[chrom] == bam_depths[chrom]
    assert _junction_signature(sam_junctions) == _junction_signature(bam_junctions)


def test_get_sam_read_data_supports_cram_with_reference(tmp_path: Path):
    """CRAM reads should load when an explicit reference FASTA is provided."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=24)

    bam_depths, bam_junctions = getSamReadData(str(fixture.bam))
    cram_depths, cram_junctions = getSamReadData(
        str(fixture.cram),
        reference_fasta=str(fixture.reference_fasta),
    )

    chrom = fixture.chrom
    assert cram_depths[chrom] == bam_depths[chrom]
    assert _junction_signature(cram_junctions) == _junction_signature(bam_junctions)


def test_get_sam_read_data_cram_with_bad_reference_is_actionable(tmp_path: Path):
    """CRAM decoding should either succeed or fail with a clear reference-related error."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)
    bad_reference = tmp_path / "missing.fa"

    try:
        depths, junctions = getSamReadData(str(fixture.cram), reference_fasta=str(bad_reference))
        assert depths
        assert fixture.chrom in depths
        assert fixture.chrom in junctions
    except Exception as exc:
        message = str(exc).lower()
        assert "reference" in message or "cram" in message or "no such file" in message


def test_get_sam_depths_and_junctions_respect_chromosome_filter(tmp_path: Path):
    """Depth and junction helpers should respect the chromosome filter argument."""
    fixture = build_alignment_fixture(tmp_path, repeat_scale=16)

    depths = getSamDepths(str(fixture.bam), chromosomes=[fixture.chrom])
    junctions = getSamJunctions(str(fixture.bam), chromosomes=[fixture.chrom], minjct=1)

    assert set(depths.keys()) == {fixture.chrom}
    assert set(junctions.keys()) == {fixture.chrom}
    assert junctions[fixture.chrom]
