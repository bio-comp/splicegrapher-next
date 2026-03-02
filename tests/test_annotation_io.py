from pathlib import Path

from SpliceGrapher.formats.annotation_io import load_gene_models
from tests.helpers.idiffir_fixture_builder import build_fixture


def test_load_gene_models_reads_gff3_and_gtf(tmp_path: Path) -> None:
    """The new annotation backend should load both GFF3 and GTF inputs."""
    fixture = build_fixture(tmp_path)
    gff3_model = load_gene_models(str(fixture.gff3))
    gtf_model = load_gene_models(str(fixture.gtf))

    assert {gene.id for gene in gff3_model.getAllGenes()} == {"GENE1", "GENE2"}
    assert {gene.id for gene in gtf_model.getAllGenes()} == {"GENE_GTF"}


def test_load_gene_models_writes_intron_cache(tmp_path: Path) -> None:
    """Passing ``outdir`` should materialize the intron BED cache artifact."""
    fixture = build_fixture(tmp_path)
    out_dir = tmp_path / "out"
    model = load_gene_models(str(fixture.gff3), outdir=out_dir)

    cache_path = out_dir / ".idiffir_cache" / "idiffir_introns.bed"
    assert model.getAllGenes()
    assert cache_path.exists()

    rows = [
        line.strip().split("\t")
        for line in cache_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert rows
    assert all(len(row) == 7 for row in rows)
    assert {"GENE1", "GENE2"} == {row[3] for row in rows}


def test_load_gene_models_normalizes_transcript_alias(tmp_path: Path) -> None:
    """`transcript` featuretype aliases should normalize to canonical mRNA."""
    gff_path = tmp_path / "alias.gff3"
    gff_path.write_text(
        "\n".join(
            [
                "chr1\ttest\tgene\t1\t100\t.\t+\t.\tID=gene1;Name=gene1",
                "chr1\ttest\ttranscript\t1\t100\t.\t+\t.\tID=tx1;Parent=gene1",
                "chr1\ttest\texon\t1\t100\t.\t+\t.\tID=ex1;Parent=tx1;gene_id=gene1;transcript_id=tx1",
                "",
            ]
        ),
        encoding="utf-8",
    )

    model = load_gene_models(str(gff_path))
    genes = model.getAllGenes()
    assert len(genes) == 1
    assert genes[0].id == "GENE1"
    assert "tx1" in genes[0].isoforms
