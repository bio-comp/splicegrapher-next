from pathlib import Path

from SpliceGrapher.formats.annotation_io import load_gene_models
from SpliceGrapher.formats.loader import loadGeneModels
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


def test_loader_wrapper_uses_annotation_backend(tmp_path: Path) -> None:
    """The legacy ``loadGeneModels`` API should delegate to annotation I/O."""
    fixture = build_fixture(tmp_path)
    out_dir = tmp_path / "wrapped_out"
    model = loadGeneModels(str(fixture.gff3), outdir=out_dir)

    assert {gene.id for gene in model.getAllGenes()} == {"GENE1", "GENE2"}
    assert (out_dir / ".idiffir_cache" / "idiffir_introns.bed").exists()
