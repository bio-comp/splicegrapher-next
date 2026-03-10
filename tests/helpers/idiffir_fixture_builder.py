from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class IDiffIRFixture:
    gff3: Path
    gtf: Path
    factor1_bam: Path
    factor2_bam: Path


def _write_models(tmp_path: Path) -> tuple[Path, Path]:
    gff3 = tmp_path / "model.gff3"
    gtf = tmp_path / "model.gtf"

    gff3.write_text(
        "##gff-version 3\n"
        "chr1\t.\tgene\t1\t300\t.\t+\t.\tID=gene1;Name=gene1\n"
        "chr1\t.\tmRNA\t1\t300\t.\t+\t.\tID=gene1.1;Parent=gene1\n"
        "chr1\t.\texon\t1\t100\t.\t+\t.\tParent=gene1.1\n"
        "chr1\t.\texon\t201\t300\t.\t+\t.\tParent=gene1.1\n"
        "chr1\t.\tgene\t401\t700\t.\t+\t.\tID=gene2;Name=gene2\n"
        "chr1\t.\tmRNA\t401\t700\t.\t+\t.\tID=gene2.1;Parent=gene2\n"
        "chr1\t.\texon\t401\t500\t.\t+\t.\tParent=gene2.1\n"
        "chr1\t.\texon\t601\t700\t.\t+\t.\tParent=gene2.1\n",
        encoding="utf-8",
    )

    gtf.write_text(
        "chr1\ttest\tgene\t1\t300\t.\t+\t.\tID=gene_gtf;Name=gene_gtf\n"
        "chr1\ttest\tmRNA\t1\t300\t.\t+\t.\tID=gene_gtf.1;Parent=gene_gtf\n"
        "chr1\ttest\texon\t1\t100\t.\t+\t.\tParent=gene_gtf.1\n"
        "chr1\ttest\texon\t201\t300\t.\t+\t.\tParent=gene_gtf.1\n",
        encoding="utf-8",
    )
    return gff3, gtf


def _segment(
    name: str, start: int, cigar: tuple[tuple[int, int], ...], query_len: int
) -> pysam.AlignedSegment:
    seg = pysam.AlignedSegment()
    seg.query_name = name
    seg.query_sequence = "A" * query_len
    seg.flag = 0
    seg.reference_id = 0
    seg.reference_start = start
    seg.mapping_quality = 60
    seg.cigartuples = list(cigar)
    seg.next_reference_id = -1
    seg.next_reference_start = -1
    seg.template_length = 0
    seg.query_qualities = pysam.qualitystring_to_array("I" * query_len)
    return seg


def _add_gene_reads(
    reads: list[tuple[int, tuple[tuple[int, int], ...], int]],
    exon1_start: int,
    intron_start: int,
    exon2_start: int,
    intron_reads: int,
    junction_reads: int,
    exon_reads: int,
) -> None:
    reads.extend([(exon1_start, ((0, 50),), 50)] * (exon_reads // 2))
    reads.extend([(exon2_start, ((0, 50),), 50)] * (exon_reads - (exon_reads // 2)))
    reads.extend([(intron_start, ((0, 50),), 50)] * intron_reads)
    reads.extend([(exon1_start + 40, ((0, 40), (3, 100), (0, 40)), 80)] * junction_reads)


def _write_bam(
    path: Path, gene1_counts: tuple[int, int, int], gene2_counts: tuple[int, int, int]
) -> None:
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": 2000}]}
    reads: list[tuple[int, tuple[tuple[int, int], ...], int]] = []
    _add_gene_reads(reads, 20, 130, 220, *gene1_counts)
    _add_gene_reads(reads, 420, 530, 620, *gene2_counts)
    reads.sort(key=lambda record: record[0])

    with pysam.AlignmentFile(str(path), "wb", header=header) as out:
        for idx, (start, cigar, query_len) in enumerate(reads):
            out.write(_segment(f"r{idx}", start, cigar, query_len))
    pysam.index(str(path))


def build_fixture(tmp_path: Path) -> IDiffIRFixture:
    gff3, gtf = _write_models(tmp_path)
    factor1_bam = tmp_path / "f1.bam"
    factor2_bam = tmp_path / "f2.bam"
    _write_bam(factor1_bam, gene1_counts=(60, 20, 20), gene2_counts=(40, 25, 20))
    _write_bam(factor2_bam, gene1_counts=(5, 40, 20), gene2_counts=(10, 35, 20))
    return IDiffIRFixture(gff3=gff3, gtf=gtf, factor1_bam=factor1_bam, factor2_bam=factor2_bam)
