from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from SpliceGrapher.formats.GeneModel import Exon, Gene, GeneModel, featureSearch


@dataclass
class _DummyGene:
    minpos: int
    maxpos: int
    strand: str = "+"

    def contains(self, loc: int, strand: str) -> bool:
        return strand == self.strand and self.minpos <= loc <= self.maxpos


def test_feature_search_supports_recursive_midpoint_indexing() -> None:
    features = [
        Exon(1, 10, "chr1", "+"),
        Exon(20, 30, "chr1", "+"),
        Exon(40, 50, "chr1", "+"),
    ]
    query = Exon(21, 22, "chr1", "+")

    assert featureSearch(features, query) is features[1]


def test_binary_search_genes_handles_large_gene_lists() -> None:
    model = GeneModel(None)
    genes = [_DummyGene(i * 10, (i * 10) + 9) for i in range(30)]

    low_gene, high_gene = model.binarySearchGenes(genes, 155, 0, len(genes) - 1)

    assert low_gene is genes[15]
    assert high_gene is genes[15]


def test_write_gff_writes_chromosome_records(tmp_path: Path) -> None:
    model = GeneModel(None)
    model.addChromosome(1, 100, "chr1")
    model.addGene(
        Gene(
            id="GENE1",
            note=None,
            start=10,
            end=20,
            chromosome="chr1",
            strand="+",
            name="GENE1",
            attr={},
        )
    )

    gff_path = tmp_path / "model.gff3"
    model.writeGFF(str(gff_path))

    lines = gff_path.read_text(encoding="utf-8").splitlines()
    expected_prefix = "Chr1\tSpliceGrapher\tchromosome\t1\t100\t.\t.\t.\tID=Chr1;Name=Chr1"
    assert lines[0].startswith(expected_prefix)
