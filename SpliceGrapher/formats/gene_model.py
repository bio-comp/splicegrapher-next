"""Gene model orchestration facade.

Domain feature entities and helper constants/functions live in
``SpliceGrapher.formats.models``. This module keeps the ``GeneModel``
container/query API and re-exports domain symbols for compatibility.
"""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import TextIO, TypeVar
from urllib.parse import unquote

from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats import models as _models
from SpliceGrapher.formats.writers.gene_model import (
    write_gff as write_gene_model_gff,
)
from SpliceGrapher.formats.writers.gene_model import (
    write_gtf as write_gene_model_gtf,
)
from SpliceGrapher.shared.progress import ProgressIndicator

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
AttrT = TypeVar("AttrT")

IntervalQuery = _models.IntervalQuery
ChromosomeGeneIndex = _models.ChromosomeGeneIndex
Chromosome = _models.Chromosome
BaseFeature = _models.BaseFeature
TranscriptRegion = _models.TranscriptRegion
Exon = _models.Exon
Isoform = _models.Isoform
CDS = _models.CDS
FpUtr = _models.FpUtr
TpUtr = _models.TpUtr
Mrna = _models.Mrna
Gene = _models.Gene
PseudoGene = _models.PseudoGene

gene_type_filter = _models.gene_type_filter
default_gene_filter = _models.default_gene_filter
dict_to_gff = _models.dict_to_gff
dict_to_gtf = _models.dict_to_gtf
cds_factory = _models.cds_factory
feature_cmp = _models.feature_cmp
feature_sort_key = _models.feature_sort_key
gene_sort_key = _models.gene_sort_key
gtf_feature_sort_key = _models.gtf_feature_sort_key
feature_overlaps = _models.feature_overlaps
feature_contains = _models.feature_contains
feature_search = _models.feature_search


@dataclass(slots=True)
class GeneModel:
    gff_path: GffRecordSource | None
    require_notes: bool = False
    chromosomes: Sequence[str] | str | None = None
    verbose: bool = False
    ignore_errors: bool = False
    # 2-dimensional model indexed by chromosome then gene
    all_genes: dict[str, Gene] = field(init=False, default_factory=dict)
    all_chr: dict[str, Chromosome] = field(init=False, default_factory=dict)
    found_types: dict[RecordType, bool] = field(init=False, default_factory=dict)
    model: dict[str, dict[str, Gene]] = field(init=False, default_factory=dict)
    mrna_forms: dict[str, dict[str, Mrna]] = field(init=False, default_factory=dict)
    mrna_gene: dict[str, dict[str, Gene]] = field(init=False, default_factory=dict)
    chromosome_index: dict[str, ChromosomeGeneIndex] = field(init=False, default_factory=dict)

    def __post_init__(self) -> None:
        """Instantiates a GeneModel object and optionally loads a GFF source."""
        if not self.gff_path:
            return
        if isinstance(self.gff_path, str) and not os.path.exists(self.gff_path):
            raise ValueError(f"Gene model file not found: {self.gff_path}")

        self.load_gene_model(
            self.gff_path,
            require_notes=self.require_notes,
            chromosomes=self.chromosomes,
            verbose=self.verbose,
            ignore_errors=self.ignore_errors,
        )

        if not self.model:
            raise ValueError(f"No gene models found in {self.gff_path}")
        self.make_sorted_model()

    def __contains__(self, gene: str | Gene) -> bool:
        """Returns true if a gene is in the model; false otherwise."""
        return str(gene) in self.all_genes

    def add_chromosome(self, start: int, end: int, name: str) -> None:
        """Adds a chromosome to a gene model or updates the end points
        if the record already exists."""
        key = name.lower()
        try:
            rec = self.all_chr[key]
            rec.minpos = min(rec.minpos, start)
            rec.maxpos = max(rec.maxpos, end)
        except KeyError:
            self.all_chr[key] = Chromosome(start, end, name)
            self.model.setdefault(key, {})

    def add_gene(self, gene: Gene) -> None:
        """Adds a gene to a gene model.  Raises a ValueError if the
        gene has already been added."""
        if gene.id in self.model[gene.chromosome] or str(gene) in self.all_genes:
            raise ValueError(f"Gene {gene.id} already stored in gene model")

        self.model[gene.chromosome][gene.id] = gene
        self.all_genes[str(gene)] = gene

    def clean_name(self, s: str) -> str:
        """
        Some feature names include URL characters that we may wish to fix.
        """
        revised = unquote(s)
        revised = revised.replace(",", "")
        return revised.replace(" ", "-")

    def get_all_acceptors(
        self, gene_filter: GeneFilter = default_gene_filter
    ) -> dict[str, dict[str, set[int]]]:
        """
        Returns a dictionary of all known acceptor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model:
            result[chrom] = self.get_known_acceptors(chrom, gene_filter)
        return result

    def get_all_donors(
        self, gene_filter: GeneFilter = default_gene_filter
    ) -> dict[str, dict[str, set[int]]]:
        """
        Returns a dictionary of all known donor sites in the gene model,
        indexed by chromosome and strand.
        """
        result = {}
        for chrom in self.model:
            result[chrom] = self.get_known_donors(chrom, gene_filter)
        return result

    def get_all_gene_ids(self, gene_filter: GeneFilter = default_gene_filter) -> list[str]:
        """Returns a list of ids for all genes stored."""
        return [g.id for g in self.all_genes.values() if gene_filter(g)]

    def get_all_genes(
        self,
        gene_filter: GeneFilter = default_gene_filter,
        *,
        verbose: bool = False,
    ) -> list[Gene]:
        """Returns a list of all genes stored."""
        indicator = ProgressIndicator(10000, verbose=verbose)
        result = []
        for g in self.all_genes.values():
            indicator.update()
            if gene_filter(g):
                result.append(g)
        indicator.finish()
        return result

    def get_annotation(
        self, key: str, annot_dict: dict[str, str], default: AttrT | None = None
    ) -> str | AttrT | None:
        """
        Convenience method for retrieving a value from an annotation dictionary
        """
        try:
            return annot_dict[key]
        except KeyError:
            return default

    def get_annotation_dict(self, s: str) -> dict[str, str]:
        """
        Parses a ';'-separated annotation string containing key-value pairs
        and returns them as a dictionary.
        """
        val_str = s.replace(" ", "")
        parts = val_str.split(";")
        result = {}
        for p in parts:
            if "=" not in p:
                continue
            key, value = p.split("=", 1)
            if not key:
                continue
            result[key] = value
        return result

    def get_chromosome(self, chr_name: str) -> Chromosome | None:
        """Returns a simple record with basic chromosome information."""
        try:
            return self.all_chr[chr_name]
        except KeyError:
            return None

    def get_chromosomes(self) -> Iterable[str]:
        """Returns a list of all chromosomes represented in the model."""
        return self.model.keys()

    def get_feature_list(self, feature_type: str | RecordType) -> list[BaseFeature]:
        """
        Returns a list of all features of the given type found in all genes.
        """
        result = []
        for gene in self.all_genes.values():
            feature_list = gene.get_feature_list(feature_type)
            if feature_list:
                result += feature_list
        return result

    def get_gene(self, chrom: str, gene_id: str) -> Gene | None:
        """
        Returns a gene from within a chromosome.  The gene will contain
        information on all exons within it.
        """
        try:
            return self.model[chrom.lower()][gene_id]
        except KeyError:
            return None

    def get_gene_by_name(self, id: str) -> Gene | None:
        """
        Returns a gene with the given id if it exists.
        """
        for k in self.model:
            try:
                return self.model[k][id.upper()]
            except KeyError:
                pass
        return None

    def get_gene_from_locations(
        self, chrom: str, start_pos: int, end_pos: int, strand: str
    ) -> Gene | None:
        """
        Finds the gene within the given chromosome that contains the given start
        and end positions.

        Uses the in-memory interval index built by ``make_sorted_model``.
        """
        chrom_key = chrom.lower()
        chrom_index = self.chromosome_index.get(chrom_key)
        if chrom_index is None:
            raise KeyError(f"Key {chrom_key} not found in {','.join(self.model.keys())}")
        return chrom_index.find_gene(start_pos, end_pos, strand)

    def get_gene_records(
        self,
        chrom: str,
        gene_filter: GeneFilter = default_gene_filter,
        *,
        verbose: bool = False,
    ) -> list[Gene]:
        """
        Returns a list of all gene instances represented within a given chromosome.
        The gene list may be filtered by changing the gene_filter function.
        """
        try:
            # return [g for g in self.model[chrom.lower()].values() if gene_filter(g)]
            indicator = ProgressIndicator(10000, verbose=verbose)
            result = []
            for g in self.model[chrom.lower()].values():
                indicator.update()
                if gene_filter(g):
                    result.append(g)
            indicator.finish()
            return result
        except KeyError:
            return []

    def get_genes(self, chrom: str) -> list[str]:
        """
        Returns a list of all genes represented within a given chromosome.
        """
        try:
            return list(self.model[chrom.lower()].keys())
        except KeyError:
            return []

    def get_genes_in_range(
        self, chrom: str, minpos: int, maxpos: int, strand: str | None = None
    ) -> list[Gene]:
        """
        Returns a list of all gene instances represented within a given chromosome
        that overlap the given range.  If no strand is specified, this will
        return all genes on both strands that overlap the range.
        """
        result = []
        for g in self.get_gene_records(chrom):
            if g.maxpos < minpos or g.minpos > maxpos:
                continue
            if strand and g.strand != strand:
                continue
            result.append(g)
        return result

    def get_known_acceptors(
        self, chrom: str, gene_filter: GeneFilter = default_gene_filter
    ) -> dict[str, set[int]]:
        """
        Returns a dictionary of all known acceptor sites for the chromosome,
        indexed by strand.
        """
        result: dict[str, set[int]] = {"-": set(), "+": set()}
        for g in self.get_gene_records(chrom, gene_filter):
            result[g.strand].update(g.acceptor_list())
        return result

    def get_known_donors(
        self, chrom: str, gene_filter: GeneFilter = default_gene_filter
    ) -> dict[str, set[int]]:
        """
        Returns a dictionary of all known donor sites for the chromosome,
        indexed by strand.
        """
        result: dict[str, set[int]] = {"-": set(), "+": set()}
        for g in self.get_gene_records(chrom, gene_filter):
            result[g.strand].update(g.donor_list())
        return result

    def get_parent(
        self,
        s: str,
        chrom: str,
        search_genes: bool = True,
        search_mrna: bool = True,
    ) -> Gene | Mrna | None:
        """
        Parent identifiers are not stored in a consistent manner.  We may have
        'AT1G01160', 'AT1G01160.1' or '12345.AT1G01160' or possibly something else.
        This method looks for the most specific candidate name to identify a record's parent.
        """
        parent_key = s.upper()
        chrom_key = chrom.lower()

        if search_mrna:
            chrom_forms = self.mrna_forms.get(chrom_key)
            if chrom_forms is not None and parent_key in chrom_forms:
                return chrom_forms[parent_key]

        if search_genes:
            chrom_genes = self.model.get(chrom_key)
            if chrom_genes is not None and parent_key in chrom_genes:
                return chrom_genes[parent_key]
        return None

    def get_record_types(self) -> list[RecordType]:
        """Returns a list of all record types found in the input file."""
        return [k for k in self.found_types if self.found_types[k]]

    def get_mrna_parent(self, chrom: str, mrna_id: str) -> Gene | None:
        """Return owning gene for a transcript identifier on ``chrom``."""
        return self.mrna_gene.get(chrom.lower(), {}).get(mrna_id.upper())

    def isoform_dict(
        self,
        gene_filter: GeneFilter = default_gene_filter,
        *,
        verbose: bool = False,
    ) -> dict[str, set[str]]:
        """Returns a dictionary that maps gene names to their corresponding
        isoform identifiers.  Each gene is associated with a set of isoform ids."""
        result: dict[str, set[str]] = {}
        for g in self.get_all_genes(gene_filter=gene_filter, verbose=verbose):
            result[g.id] = set(list(g.mrna.keys()) + list(g.isoforms.keys()))
        return result

    def load_gene_model(
        self,
        gff_records: GffRecordSource,
        *,
        require_notes: bool = False,
        chromosomes: Sequence[str] | str | None = None,
        verbose: bool = False,
        ignore_errors: bool = False,
    ) -> None:
        """
        Reads a tab-delimited gene annotation GFF file and stores information
        on chromosomes, the genes within each chromosome and exons within each gene.

        Parameters:
          'gff_records'  - source of GFF records; may be a file path, a file stream
                           or a list/set of strings
          'require_notes' - require annotations for all gene records (default=False)
          'chromosomes'  - chromosome name or list of chromosomes to store (default=all)
          'verbose'      - provide verbose feedback (default=False)
          'ignore_errors' - ignore error conditions (default=False)
        """
        from SpliceGrapher.formats.parsers.gene_model_gff import load_gene_model_records

        load_gene_model_records(
            self,
            gff_records,
            require_notes=require_notes,
            chromosomes=chromosomes,
            verbose=verbose,
            ignore_errors=ignore_errors,
        )

    def make_sorted_model(self) -> None:
        self.chromosome_index = {
            chrom: ChromosomeGeneIndex.build(self.model[chrom].values()) for chrom in self.model
        }

    def write_gff(
        self,
        gff_path: str | TextIO,
        *,
        gene_filter: GeneFilter = default_gene_filter,
        gene_set: set[str] | list[str] | tuple[str, ...] | None = None,
        verbose: bool = False,
    ) -> None:
        """Writes a complete gene model out to a GFF file."""
        write_gene_model_gff(
            self,
            gff_path,
            gene_filter=gene_filter,
            gene_set=gene_set,
            verbose=verbose,
        )

    def write_gtf(
        self,
        gtf_path: str | TextIO,
        *,
        gene_filter: GeneFilter = default_gene_filter,
        verbose: bool = False,
    ) -> None:
        """Writes a complete gene model out to a GTF file."""
        write_gene_model_gtf(
            self,
            gtf_path,
            gene_filter=gene_filter,
            verbose=verbose,
        )
