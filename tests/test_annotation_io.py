from pathlib import Path
from typing import Callable, cast

from SpliceGrapher.formats.annotation_io import load_gene_models
from tests.helpers.idiffir_fixture_builder import build_fixture


def test_load_gene_models_reads_gff3_and_gtf(tmp_path: Path) -> None:
    """The new annotation backend should load both GFF3 and GTF inputs."""
    fixture = build_fixture(tmp_path)
    gff3_model = load_gene_models(str(fixture.gff3))
    gtf_model = load_gene_models(str(fixture.gtf))

    assert {gene.id for gene in gff3_model.get_all_genes()} == {"GENE1", "GENE2"}
    assert {gene.id for gene in gtf_model.get_all_genes()} == {"GENE_GTF"}


def test_load_gene_models_writes_intron_cache(tmp_path: Path) -> None:
    """Passing ``outdir`` should materialize the intron BED cache artifact."""
    fixture = build_fixture(tmp_path)
    out_dir = tmp_path / "out"
    model = load_gene_models(str(fixture.gff3), outdir=out_dir)

    cache_path = out_dir / ".idiffir_cache" / "idiffir_introns.bed"
    assert model.get_all_genes()
    assert cache_path.exists()

    rows = [
        line.strip().split("\t")
        for line in cache_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    assert rows
    assert all(len(row) == 7 for row in rows)
    assert {"GENE1", "GENE2"} == {row[3] for row in rows}


def test_load_gene_models_writes_gffutils_cache_db(tmp_path: Path) -> None:
    """Passing ``cache_dir`` should materialize the deterministic sqlite cache."""
    fixture = build_fixture(tmp_path)
    cache_dir = tmp_path / "db-cache"

    model = load_gene_models(str(fixture.gff3), cache_dir=cache_dir)

    assert model.get_all_genes()
    db_files = sorted(cache_dir.glob("annotation_*.db"))
    assert len(db_files) == 1


def test_load_gene_models_rejects_unknown_keyword_argument(tmp_path: Path) -> None:
    """The loader boundary should reject stray keyword arguments immediately."""
    fixture = build_fixture(tmp_path)
    call_loader = cast(Callable[..., object], load_gene_models)

    try:
        call_loader(str(fixture.gff3), nonsense=True)
    except TypeError as exc:
        assert "nonsense" in str(exc)
    else:
        raise AssertionError("Expected TypeError for unsupported keyword argument")


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
    genes = model.get_all_genes()
    assert len(genes) == 1
    assert genes[0].id == "GENE1"
    assert "tx1" in genes[0].transcripts
