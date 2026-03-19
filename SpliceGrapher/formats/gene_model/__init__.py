"""Gene model orchestration facade.

Domain feature entities and helper constants/functions live in
`SpliceGrapher.formats.models`. This package keeps the `GeneModel`
container/query API and re-exports domain symbols for compatibility.
"""

from __future__ import annotations

from SpliceGrapher.formats import models as _models
from SpliceGrapher.formats.writers.gene_model import write_gff as write_gene_model_gff
from SpliceGrapher.formats.writers.gene_model import write_gtf as write_gene_model_gtf

from .model import GeneModel

KNOWN_RECTYPES = _models.KNOWN_RECTYPES
IGNORE_RECTYPES = _models.IGNORE_RECTYPES
CDS_TYPES = _models.CDS_TYPES
GTF_GENE_ID = _models.GTF_GENE_ID
GTF_GENE_NAME = _models.GTF_GENE_NAME
GTF_TRANSCRIPT = _models.GTF_TRANSCRIPT
GTF_TRANSNAME = _models.GTF_TRANSNAME
GTF_SOURCE = _models.GTF_SOURCE
GTF_EXON_ID = _models.GTF_EXON_ID
GTF_PROTEIN_ID = _models.GTF_PROTEIN_ID
RECTYPE_MAP = _models.RECTYPE_MAP
ISOFORM_TYPE = _models.ISOFORM_TYPE
ID_FIELD = _models.ID_FIELD
NAME_FIELD = _models.NAME_FIELD
NOTE_FIELD = _models.NOTE_FIELD
PARENT_FIELD = _models.PARENT_FIELD
POSSIBLE_GENE_FIELDS = _models.POSSIBLE_GENE_FIELDS
POSSIBLE_FORM_FIELDS = _models.POSSIBLE_FORM_FIELDS
GFF_ID = _models.GFF_ID
MAX_BAD_LINES = _models.MAX_BAD_LINES
VALID_STRANDS = _models.VALID_STRANDS
FORM_DELIMITERS = _models.FORM_DELIMITERS
SPLICE_DIMER_OFFSET = _models.SPLICE_DIMER_OFFSET
GTF_ORDER_POLICY = _models.GTF_ORDER_POLICY

GeneFilter = _models.GeneFilter
GffRecordSource = _models.GffRecordSource
AttrValue = _models.AttrValue
AttrMap = _models.AttrMap
ExtraAttrMap = _models.ExtraAttrMap

IntervalQuery = _models.IntervalQuery
ChromosomeGeneIndex = _models.ChromosomeGeneIndex
Chromosome = _models.Chromosome
BaseFeature = _models.BaseFeature
TranscriptRegion = _models.TranscriptRegion
Exon = _models.Exon
Transcript = _models.Transcript
CDS = _models.CDS
FpUtr = _models.FpUtr
TpUtr = _models.TpUtr
Gene = _models.Gene
PseudoGene = _models.PseudoGene

gene_type_filter = _models.gene_type_filter
default_gene_filter = _models.default_gene_filter
cds_factory = _models.cds_factory
feature_cmp = _models.feature_cmp
feature_sort_key = _models.feature_sort_key
gene_sort_key = _models.gene_sort_key
gtf_feature_sort_key = _models.gtf_feature_sort_key
feature_overlaps = _models.feature_overlaps
feature_contains = _models.feature_contains

__all__ = [
    "AttrMap",
    "AttrValue",
    "BaseFeature",
    "CDS",
    "CDS_TYPES",
    "Chromosome",
    "ChromosomeGeneIndex",
    "Exon",
    "ExtraAttrMap",
    "FORM_DELIMITERS",
    "FpUtr",
    "GFF_ID",
    "GTF_EXON_ID",
    "GTF_GENE_ID",
    "GTF_GENE_NAME",
    "GTF_ORDER_POLICY",
    "GTF_PROTEIN_ID",
    "GTF_SOURCE",
    "GTF_TRANSCRIPT",
    "GTF_TRANSNAME",
    "Gene",
    "GeneFilter",
    "GeneModel",
    "GffRecordSource",
    "ID_FIELD",
    "IGNORE_RECTYPES",
    "ISOFORM_TYPE",
    "IntervalQuery",
    "KNOWN_RECTYPES",
    "MAX_BAD_LINES",
    "NAME_FIELD",
    "NOTE_FIELD",
    "PARENT_FIELD",
    "POSSIBLE_FORM_FIELDS",
    "POSSIBLE_GENE_FIELDS",
    "PseudoGene",
    "RECTYPE_MAP",
    "SPLICE_DIMER_OFFSET",
    "TpUtr",
    "Transcript",
    "TranscriptRegion",
    "VALID_STRANDS",
    "cds_factory",
    "default_gene_filter",
    "feature_cmp",
    "feature_contains",
    "feature_overlaps",
    "feature_sort_key",
    "gene_sort_key",
    "gene_type_filter",
    "gtf_feature_sort_key",
    "write_gene_model_gff",
    "write_gene_model_gtf",
]
