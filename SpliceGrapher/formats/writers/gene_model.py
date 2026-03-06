"""Writer boundary for ``GeneModel`` serialization."""

from __future__ import annotations

from collections.abc import Callable, Collection, Iterator
from contextlib import contextmanager
from typing import TYPE_CHECKING, TextIO

from SpliceGrapher.formats.serializers import (
    gff_line_for_chromosome,
    gff_text_for_gene,
    gtf_text_for_gene,
)
from SpliceGrapher.shared.progress import ProgressIndicator

if TYPE_CHECKING:
    from SpliceGrapher.formats.gene_model import Gene, GeneModel

GeneFilter = Callable[["Gene"], bool]


def _default_gene_filter(gene: Gene) -> bool:
    return True


def _gene_sort_key(gene: Gene) -> tuple[int, int, str]:
    return (gene.minpos, gene.maxpos, gene.id)


@contextmanager
def _open_output(path: str | TextIO) -> Iterator[TextIO]:
    if isinstance(path, str):
        with open(path, "w") as stream:
            yield stream
        return
    yield path


def write_gff(
    model: GeneModel,
    gff_path: str | TextIO,
    *,
    gene_filter: GeneFilter = _default_gene_filter,
    gene_set: Collection[str] | None = None,
    verbose: bool = False,
) -> None:
    """Write full gene model to GFF output."""
    with _open_output(gff_path) as out_stream:
        chrom_list = sorted(model.all_chr)
        indicator = ProgressIndicator(10000, verbose=verbose)
        for chrom_name in chrom_list:
            chrom = model.get_chromosome(chrom_name)
            if chrom is None:
                continue
            out_stream.write(f"{gff_line_for_chromosome(chrom)}\n")
            genes = model.get_gene_records(chrom_name, gene_filter)
            if gene_set:
                genes = [g for g in genes if g.id in gene_set or g.name in gene_set]
            genes.sort(key=_gene_sort_key)
            for gene in genes:
                indicator.update()
                strings = gff_text_for_gene(gene)
                if strings:
                    out_stream.write(f"{strings}\n")
        indicator.finish()


def write_gtf(
    model: GeneModel,
    gtf_path: str | TextIO,
    *,
    gene_filter: GeneFilter = _default_gene_filter,
    verbose: bool = False,
) -> None:
    """Write full gene model to GTF output."""
    with _open_output(gtf_path) as out_stream:
        chrom_list = sorted(model.all_chr)
        indicator = ProgressIndicator(10000, verbose=verbose)
        for chrom_name in chrom_list:
            genes = model.get_gene_records(chrom_name, gene_filter)
            genes.sort(key=_gene_sort_key)
            for gene in genes:
                indicator.update()
                strings = gtf_text_for_gene(gene)
                if strings:
                    out_stream.write(f"{strings}\n")
        indicator.finish()


__all__ = ["write_gff", "write_gtf"]
