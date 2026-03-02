"""Writer boundary for ``GeneModel`` serialization."""

from __future__ import annotations

from typing import TYPE_CHECKING, TextIO

from SpliceGrapher.formats.GeneModel import defaultGeneFilter, geneSortKey
from SpliceGrapher.shared.process_utils import getAttribute
from SpliceGrapher.shared.progress import ProgressIndicator

if TYPE_CHECKING:
    from SpliceGrapher.formats.GeneModel import GeneFilter, GeneModel


def write_gff(model: "GeneModel", gffPath: str | TextIO, **args: object) -> None:
    """Write full gene model to GFF output."""
    geneFilter = getAttribute("geneFilter", defaultGeneFilter, **args)
    geneSubset = getAttribute("geneSet", None, **args)
    verbose = getAttribute("verbose", False, **args)

    outStream = open(gffPath, "w") if isinstance(gffPath, str) else gffPath
    chromList = sorted(model.allChr.keys())
    indicator = ProgressIndicator(10000, verbose=verbose)
    for c in chromList:
        chrom = model.get_chromosome(c)
        if chrom is None:
            continue
        outStream.write(f"{chrom.gffString()}\n")
        genes = model.get_gene_records(c, geneFilter)
        if geneSubset:
            genes = [g for g in genes if g.id in geneSubset or g.name in geneSubset]
        genes.sort(key=geneSortKey)
        for g in genes:
            indicator.update()
            strings = g.gffStrings()
            if strings:
                outStream.write(f"{strings}\n")
    indicator.finish()


def write_gtf(
    model: "GeneModel",
    gtfPath: str | TextIO,
    geneFilter: "GeneFilter" = defaultGeneFilter,
    **args: object,
) -> None:
    """Write full gene model to GTF output."""
    verbose = getAttribute("verbose", False, **args)
    outStream = open(gtfPath, "w") if isinstance(gtfPath, str) else gtfPath
    chromList = sorted(model.allChr.keys())
    indicator = ProgressIndicator(10000, verbose=verbose)
    for c in chromList:
        genes = model.get_gene_records(c, geneFilter)
        genes.sort(key=geneSortKey)
        for g in genes:
            indicator.update()
            strings = g.gtfStrings()
            if strings:
                outStream.write(f"{strings}\n")
    indicator.finish()


__all__ = ["write_gff", "write_gtf"]
