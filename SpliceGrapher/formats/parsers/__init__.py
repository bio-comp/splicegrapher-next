"""Parser boundaries for format modules."""

from .gene_model_gff import load_gene_model_records
from .splice_graph import SpliceGraphParser

__all__ = ["SpliceGraphParser", "load_gene_model_records"]
