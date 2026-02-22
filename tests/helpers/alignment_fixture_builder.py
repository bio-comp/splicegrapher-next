from dataclasses import dataclass
from pathlib import Path

import pysam


@dataclass(frozen=True)
class AlignmentFixture:
    reference_fasta: Path
    gff3: Path
    gtf: Path
    bam: Path
    cram: Path
    sam: Path
    chrom: str
    region_start: int
    region_end: int


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
    seg.cigar = cigar
    seg.next_reference_id = -1
    seg.next_reference_start = -1
    seg.template_length = 0
    seg.query_qualities = pysam.qualitystring_to_array("I" * query_len)
    return seg


def _write_models(tmp_path: Path, chrom: str) -> tuple[Path, Path]:
    gff3 = tmp_path / "model.gff3"
    gtf = tmp_path / "model.gtf"

    gff3.write_text(
        "##gff-version 3\n"
        f"{chrom}\t.\tgene\t1\t1000\t.\t+\t.\tID=geneA;Name=geneA\n"
        f"{chrom}\t.\tmRNA\t1\t1000\t.\t+\t.\tID=geneA.1;Parent=geneA\n"
        f"{chrom}\t.\texon\t1\t300\t.\t+\t.\tParent=geneA.1\n"
        f"{chrom}\t.\texon\t601\t1000\t.\t+\t.\tParent=geneA.1\n",
        encoding="utf-8",
    )

    gtf.write_text(
        f'{chrom}\ttest\texon\t1\t300\t.\t+\t.\tgene_id "geneA"; transcript_id "geneA.1";\n'
        f'{chrom}\ttest\texon\t601\t1000\t.\t+\t.\tgene_id "geneA"; transcript_id "geneA.1";\n',
        encoding="utf-8",
    )
    return gff3, gtf


def build_alignment_fixture(tmp_path: Path, *, repeat_scale: int = 120) -> AlignmentFixture:
    chrom = "chr1"
    reference_fasta = tmp_path / "ref.fa"
    reference_fasta.write_text(f">{chrom}\n" + ("A" * 4000) + "\n", encoding="utf-8")
    pysam.faidx(str(reference_fasta))

    gff3, gtf = _write_models(tmp_path, chrom)
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": chrom, "LN": 4000}]}

    reads: list[tuple[int, tuple[tuple[int, int], ...], int]] = []
    for i in range(repeat_scale):
        # Fully exonic read.
        reads.append((40 + (i % 30), ((0, 50),), 50))
        # Spliced read with one junction.
        reads.append((120 + (i % 20), ((0, 30), (3, 100), (0, 30)), 60))
        # Read with deletion.
        reads.append((260 + (i % 15), ((0, 25), (2, 5), (0, 20)), 45))

    reads.sort(key=lambda item: item[0])
    segments = [
        _segment(f"r{idx}", start, cigar, query_len)
        for idx, (start, cigar, query_len) in enumerate(reads)
    ]

    bam = tmp_path / "reads.bam"
    with pysam.AlignmentFile(bam, "wb", header=header) as out:
        for seg in segments:
            out.write(seg)
    pysam.index(str(bam))

    sam = tmp_path / "reads.sam"
    with pysam.AlignmentFile(sam, "w", header=header) as out:
        for seg in segments:
            out.write(seg)

    cram = tmp_path / "reads.cram"
    with pysam.AlignmentFile(
        cram, "wc", header=header, reference_filename=str(reference_fasta)
    ) as out:
        for seg in segments:
            out.write(seg)
    pysam.index(str(cram))

    return AlignmentFixture(
        reference_fasta=reference_fasta,
        gff3=gff3,
        gtf=gtf,
        bam=bam,
        cram=cram,
        sam=sam,
        chrom=chrom,
        region_start=1,
        region_end=1000,
    )
