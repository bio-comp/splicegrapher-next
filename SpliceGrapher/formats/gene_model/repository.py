"""Repository boundary for GeneModel loading and serialization."""

from __future__ import annotations

from collections.abc import Sequence
from typing import TYPE_CHECKING, TextIO

from SpliceGrapher.formats import models as _models

if TYPE_CHECKING:
    from .model import GeneModel

GeneFilter = _models.GeneFilter
GffRecordSource = _models.GffRecordSource
default_gene_filter = _models.default_gene_filter


class GeneModelRepository:
    """Serialization boundary for loading and writing ``GeneModel`` payloads."""

    @staticmethod
    def load(
        model: GeneModel,
        gff_records: GffRecordSource,
        *,
        require_notes: bool = False,
        chromosomes: Sequence[str] | str | None = None,
        verbose: bool = False,
        ignore_errors: bool = False,
    ) -> None:
        from SpliceGrapher.formats.parsers.gene_model_gff import load_gene_model_records

        load_gene_model_records(
            model,
            gff_records,
            require_notes=require_notes,
            chromosomes=chromosomes,
            verbose=verbose,
            ignore_errors=ignore_errors,
        )

    @staticmethod
    def write_gff(
        model: GeneModel,
        gff_path: str | TextIO,
        *,
        gene_filter: GeneFilter = default_gene_filter,
        gene_set: set[str] | list[str] | tuple[str, ...] | None = None,
        verbose: bool = False,
    ) -> None:
        from SpliceGrapher.formats import gene_model as gene_model_package

        gene_model_package.write_gene_model_gff(
            model,
            gff_path,
            gene_filter=gene_filter,
            gene_set=gene_set,
            verbose=verbose,
        )

    @staticmethod
    def write_gtf(
        model: GeneModel,
        gtf_path: str | TextIO,
        *,
        gene_filter: GeneFilter = default_gene_filter,
        verbose: bool = False,
    ) -> None:
        from SpliceGrapher.formats import gene_model as gene_model_package

        gene_model_package.write_gene_model_gtf(
            model,
            gtf_path,
            gene_filter=gene_filter,
            verbose=verbose,
        )
