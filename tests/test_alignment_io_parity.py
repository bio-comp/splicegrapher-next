from pathlib import Path

import numpy

from SpliceGrapher.formats.alignment_io import getSamJunctions, getSamReadData
from tests.helpers.alignment_fixture_builder import build_alignment_fixture
from tests.helpers.legacy_depth_reference import compute_depths_and_junctions


def _depth_window(
    depths_by_chrom: dict[str, list[int]], chrom: str, start: int, end: int
) -> numpy.ndarray:
    return numpy.array(depths_by_chrom[chrom][start : end + 1], dtype=int)


def _trim_expected_depths(expected_depths: numpy.ndarray, observed_length: int) -> numpy.ndarray:
    return expected_depths[:observed_length]


def _junction_dict(junctions_by_chrom, chrom: str, start: int) -> dict[tuple[int, int], int]:
    min_idx = start - 1
    result: dict[tuple[int, int], int] = {}
    for jct in junctions_by_chrom.get(chrom, []):
        donor = jct.donor() - min_idx
        acceptor = (jct.acceptor() - 1) - min_idx
        result[(donor, acceptor)] = jct.count
    return result


def test_alignment_fixture_builds_sam_bam_cram_and_annotations(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)
    assert fixture.gff3.exists()
    assert fixture.gtf.exists()
    assert fixture.bam.exists()
    assert (fixture.bam.parent / (fixture.bam.name + ".bai")).exists()
    assert fixture.sam.exists()
    assert fixture.cram.exists()
    assert (fixture.cram.parent / (fixture.cram.name + ".crai")).exists()
    assert fixture.reference_fasta.exists()


def test_bam_depth_and_junction_parity_against_legacy_reference(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    actual_depths, actual_junctions = getSamReadData(str(fixture.bam))
    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.bam, fixture.chrom, fixture.region_start, fixture.region_end
    )

    actual_depth_array = _depth_window(
        actual_depths, fixture.chrom, fixture.region_start, fixture.region_end
    )
    actual_junction_dict = _junction_dict(actual_junctions, fixture.chrom, fixture.region_start)
    trimmed_expected_depths = _trim_expected_depths(expected_depths, len(actual_depth_array))

    assert numpy.array_equal(actual_depth_array, trimmed_expected_depths)
    assert actual_junction_dict == expected_junctions


def test_cram_reference_depths_available_for_parity_harness(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    cram_depths, cram_junctions = getSamReadData(
        str(fixture.cram),
        reference_fasta=str(fixture.reference_fasta),
    )
    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.cram,
        fixture.chrom,
        fixture.region_start,
        fixture.region_end,
        reference_fasta=fixture.reference_fasta,
    )

    actual_depth_array = _depth_window(
        cram_depths, fixture.chrom, fixture.region_start, fixture.region_end
    )
    actual_junction_dict = _junction_dict(cram_junctions, fixture.chrom, fixture.region_start)
    trimmed_expected_depths = _trim_expected_depths(expected_depths, len(actual_depth_array))

    assert len(actual_depth_array) <= fixture.region_end - fixture.region_start + 1
    assert numpy.array_equal(actual_depth_array, trimmed_expected_depths)
    assert actual_junction_dict == expected_junctions


def test_sam_reference_depths_available_for_parity_harness(tmp_path: Path):
    fixture = build_alignment_fixture(tmp_path)

    sam_junctions = getSamJunctions(str(fixture.sam), minjct=1)
    sam_depths, _ = getSamReadData(str(fixture.sam), junctions=False)
    expected_depths, expected_junctions = compute_depths_and_junctions(
        fixture.sam,
        fixture.chrom,
        fixture.region_start,
        fixture.region_end,
    )

    actual_depth_array = _depth_window(
        sam_depths, fixture.chrom, fixture.region_start, fixture.region_end
    )
    actual_junction_dict = _junction_dict(sam_junctions, fixture.chrom, fixture.region_start)
    trimmed_expected_depths = _trim_expected_depths(expected_depths, len(actual_depth_array))

    assert numpy.array_equal(actual_depth_array, trimmed_expected_depths)
    assert actual_junction_dict == expected_junctions
