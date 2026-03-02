"""Writer boundary for ``GeneModel`` serialization."""

from __future__ import annotations

from collections.abc import Callable, Collection
from typing import TYPE_CHECKING, TextIO

from SpliceGrapher.shared.progress import ProgressIndicator

if TYPE_CHECKING:
    from SpliceGrapher.formats.GeneModel import Gene, GeneModel

GeneFilter = Callable[["Gene"], bool]


def _default_gene_filter(gene: "Gene") -> bool:
    return True


def _gene_sort_key(gene: "Gene") -> tuple[int, int, str]:
    return (gene.minpos, gene.maxpos, gene.id)


def write_gff(
    model: "GeneModel",
    gff_path: str | TextIO,
    *,
    gene_filter: GeneFilter = _default_gene_filter,
    gene_set: Collection[str] | None = None,
    verbose: bool = False,
) -> None:
    """Write full gene model to GFF output."""
    outStream = open(gff_path, "w") if isinstance(gff_path, str) else gff_path
    chromList = sorted(model.allChr.keys())
    indicator = ProgressIndicator(10000, verbose=verbose)
    for c in chromList:
        chrom = model.get_chromosome(c)
        if chrom is None:
            continue
        outStream.write(f"{chrom.gffString()}\n")
        genes = model.get_gene_records(c, gene_filter)
        if gene_set:
            genes = [g for g in genes if g.id in gene_set or g.name in gene_set]
        genes.sort(key=_gene_sort_key)
        for g in genes:
            indicator.update()
            strings = g.gffStrings()
            if strings:
                outStream.write(f"{strings}\n")
    indicator.finish()


def write_gtf(
    model: "GeneModel",
    gtf_path: str | TextIO,
    *,
    gene_filter: GeneFilter = _default_gene_filter,
    verbose: bool = False,
) -> None:
    """Write full gene model to GTF output."""
    outStream = open(gtf_path, "w") if isinstance(gtf_path, str) else gtf_path
    chromList = sorted(model.allChr.keys())
    indicator = ProgressIndicator(10000, verbose=verbose)
    for c in chromList:
        genes = model.get_gene_records(c, gene_filter)
        genes.sort(key=_gene_sort_key)
        for g in genes:
            indicator.update()
            strings = g.gtfStrings()
            if strings:
                outStream.write(f"{strings}\n")
    indicator.finish()


__all__ = ["write_gff", "write_gtf"]
