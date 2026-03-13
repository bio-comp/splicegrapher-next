"""GeneModel container and query API."""

from __future__ import annotations

import os
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from typing import TextIO
from urllib.parse import unquote

from SpliceGrapher.core.enums import RecordType
from SpliceGrapher.formats import models as _models

from .repository import GeneModelRepository

GeneFilter = _models.GeneFilter
GffRecordSource = _models.GffRecordSource
ChromosomeGeneIndex = _models.ChromosomeGeneIndex
Chromosome = _models.Chromosome
BaseFeature = _models.BaseFeature
Transcript = _models.Transcript
Gene = _models.Gene
default_gene_filter = _models.default_gene_filter


@dataclass(slots=True, kw_only=True)
class GeneModel:
    require_notes: bool = False
    chromosomes: Sequence[str] | str | None = None
    verbose: bool = False
    ignore_errors: bool = False
    # 2-dimensional model indexed by chromosome then gene
    all_genes: dict[str, Gene] = field(init=False, default_factory=dict)
    all_chr: dict[str, Chromosome] = field(init=False, default_factory=dict)
    found_types: dict[RecordType, bool] = field(init=False, default_factory=dict)
    model: dict[str, dict[str, Gene]] = field(init=False, default_factory=dict)
    mrna_forms: dict[str, dict[str, Transcript]] = field(init=False, default_factory=dict)
    mrna_gene: dict[str, dict[str, Gene]] = field(init=False, default_factory=dict)
    chromosome_index: dict[str, ChromosomeGeneIndex] = field(init=False, default_factory=dict)

    @classmethod
    def from_gff(
        cls,
        gff_records: GffRecordSource,
        *,
        require_notes: bool = False,
        chromosomes: Sequence[str] | str | None = None,
        verbose: bool = False,
        ignore_errors: bool = False,
    ) -> GeneModel:
        """Factory method that loads model state from GFF records."""
        if isinstance(gff_records, str) and not os.path.exists(gff_records):
            raise ValueError(f"Gene model file not found: {gff_records}")

        model = cls(
            require_notes=require_notes,
            chromosomes=chromosomes,
            verbose=verbose,
            ignore_errors=ignore_errors,
        )
        model.load_gene_model(
            gff_records,
            require_notes=require_notes,
            chromosomes=chromosomes,
            verbose=verbose,
            ignore_errors=ignore_errors,
        )
        if not model.model:
            raise ValueError(f"No gene models found in {gff_records}")
        model.make_sorted_model()
        return model

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
        chrom_key = gene.chromosome.lower()
        chrom_genes = self.model.setdefault(chrom_key, {})
        if self._find_gene_case_insensitive(chrom_genes, gene.id) or str(gene) in self.all_genes:
            raise ValueError(f"Gene {gene.id} already stored in gene model")

        self.model[chrom_key][gene.id] = gene
        self.all_genes[str(gene)] = gene

    @staticmethod
    def _find_gene_case_insensitive(chrom_genes: dict[str, Gene], gene_id: str) -> Gene | None:
        gene = chrom_genes.get(gene_id)
        if gene is not None:
            return gene
        search_id = gene_id.upper()
        gene = chrom_genes.get(search_id)
        if gene is not None:
            return gene
        for stored_id, stored_gene in chrom_genes.items():
            if stored_id.upper() == search_id:
                return stored_gene
        return None

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
        _ = verbose  # Compatibility: retained keyword parameter.
        return [g for g in self.all_genes.values() if gene_filter(g)]

    def iter_all_genes(self) -> Iterable[Gene]:
        """Yield all genes."""
        yield from self.all_genes.values()

    def get_annotation_dict(self, s: str) -> dict[str, str]:
        """
        Parses a ';'-separated annotation string containing key-value pairs
        and returns them as a dictionary.
        """
        return {
            key: value
            for part in s.replace(" ", "").split(";")
            if "=" in part
            for key, value in [part.split("=", 1)]
            if key
        }

    def get_chromosome(self, chr_name: str) -> Chromosome | None:
        """Returns a simple record with basic chromosome information."""
        return self.all_chr.get(chr_name.lower())

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
                result.extend(feature_list)
        return result

    def get_gene(self, chrom: str, gene_id: str) -> Gene | None:
        """
        Returns a gene from within a chromosome.  The gene will contain
        information on all exons within it.
        """
        chrom_genes = self.model.get(chrom.lower())
        if chrom_genes is None:
            return None
        return self._find_gene_case_insensitive(chrom_genes, gene_id)

    def get_gene_by_name(self, gene_id: str) -> Gene | None:
        """
        Returns a gene with the given id if it exists.
        """
        for chrom_genes in self.model.values():
            gene = self._find_gene_case_insensitive(chrom_genes, gene_id)
            if gene is not None:
                return gene
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
        _ = verbose  # Compatibility: retained keyword parameter.
        chrom_genes = self.model.get(chrom.lower())
        if chrom_genes is None:
            return []
        return [g for g in chrom_genes.values() if gene_filter(g)]

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
        chrom_key = chrom.lower()
        chrom_index = self.chromosome_index.get(chrom_key)
        if chrom_index is not None:
            if strand is not None:
                return chrom_index.find_genes_overlapping(minpos, maxpos, strand)
            result: list[Gene] = []
            for strand_key in ("+", "-", "."):
                result.extend(chrom_index.find_genes_overlapping(minpos, maxpos, strand_key))
            return result

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
    ) -> Gene | Transcript | None:
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
            result[g.id] = set(g.transcripts.keys())
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
        GeneModelRepository.load(
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
        GeneModelRepository.write_gff(
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
        GeneModelRepository.write_gtf(
            self,
            gtf_path,
            gene_filter=gene_filter,
            verbose=verbose,
        )
