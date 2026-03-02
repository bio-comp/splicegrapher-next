"""GFF parser boundary for ``GeneModel`` loading."""

from __future__ import annotations

import sys
from collections.abc import Iterable
from typing import TYPE_CHECKING

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import RecordType, Strand
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.format_utils import comma_format
from SpliceGrapher.shared.process_utils import getAttribute
from SpliceGrapher.shared.progress import ProgressIndicator

if TYPE_CHECKING:
    from SpliceGrapher.formats.GeneModel import GeneModel, GffRecordSource


def load_gene_model_records(
    model: "GeneModel", gffRecords: "GffRecordSource", **args: object
) -> None:
    """Load gene model records from GFF-like input into ``model``."""
    import SpliceGrapher.formats.GeneModel as gm

    verbose = getAttribute("verbose", False, **args)
    requireNotes = getAttribute("requireNodes", False, **args)
    ignoreErrors = getAttribute("ignoreErrors", False, **args)
    chromosomes: list[str] | None = None

    # Convenience method for handling exceptions that could be ignored:
    def conditionalException(message: str) -> None:
        if not ignoreErrors:
            raise RuntimeError(message)

    if isinstance(gffRecords, str):
        if verbose:
            sys.stderr.write(f"Loading and validating gene models in {gffRecords}\n")
        instream: Iterable[str] = ez_open(gffRecords)
    elif isinstance(gffRecords, (list, set, tuple)):
        if verbose:
            sys.stderr.write(f"Loading and validating gene models in {len(gffRecords)} records\n")
        instream = gffRecords
    else:
        raise ValueError(
            "Unrecognized GFF record source "
            f"({type(gffRecords).__name__}); must be file path or a list/set of strings."
        )

    if "chromosomes" in args:
        chrParam = args["chromosomes"]
        if isinstance(chrParam, (list, tuple, set)):
            chromosomes = [str(s).lower() for s in chrParam]
        elif isinstance(chrParam, str) and chrParam:
            chromosomes = [chrParam.lower()]
        if verbose and chromosomes is not None:
            sys.stderr.write(f"GeneModel loading chromosomes {','.join(chromosomes)}\n")

    model.model = {}
    model.mRNAforms = {}
    model.allGenes = {}
    model.foundTypes = {}
    model.allChr = {}
    lineCtr = 0

    geneAlias: dict[str, str] = {}
    geneCount = 0
    exonCount = 0
    isoCount = 0
    mrnaCount = 0
    cdsCount = 0
    badLines = 0
    indicator = ProgressIndicator(1000000, verbose=verbose)

    for line in instream:
        lineCtr += 1
        indicator.update()
        s = line.rstrip()
        if not s or s[0] == "#":
            continue

        parts = s.split("\t")
        if len(parts) < 7:
            badLines += 1
            if verbose:
                sys.stderr.write(
                    f"line {lineCtr}: invalid GFF format "
                    "(not enough columns); file may be corrupt\n"
                )
            if badLines >= gm.MAX_BAD_LINES:
                sys.stderr.write("\nInput GFF file appears to be invalid; aborting.\n")
                raise ValueError("Invalid GFF input file")
            continue

        annots = model.get_annotation_dict(parts[-1])
        chrName = parts[0].lower()
        if chromosomes and chrName not in chromosomes:
            continue

        try:
            recType = coerce_enum(parts[2].lower(), RecordType, field="record_type")
        except ValueError as exc:
            raise ValueError(f"line {lineCtr}: unknown record type '{parts[2]}'") from exc
        recType = gm.RECTYPE_MAP.get(recType, recType)

        startPos = int(parts[3])
        endPos = int(parts[4])
        try:
            strand = coerce_enum(parts[6], Strand, field="strand")
        except ValueError as exc:
            raise ValueError(f"line {lineCtr}: unknown strand '{parts[6]}'") from exc

        model.foundTypes[recType] = (
            recType in gm.KNOWN_RECTYPES and recType not in gm.IGNORE_RECTYPES
        )
        if not model.foundTypes[recType]:
            continue

        if recType in [RecordType.GENE, RecordType.PSEUDOGENE]:
            try:
                gid = annots[gm.ID_FIELD].upper()
            except KeyError:
                raise ValueError(f"line {lineCtr}: {recType} record has no ID field:\n{line}\n")

            name = model.get_annotation(gm.NAME_FIELD, annots, None)
            if name:
                name = model.clean_name(name)

            note = model.get_annotation(gm.NOTE_FIELD, annots)
            if not note and requireNotes:
                continue

            if strand not in gm.VALID_STRANDS:
                conditionalException(f"line {lineCtr}: {recType} record with unknown strand")

            if chrName not in model.model:
                model.model[chrName] = {}
                model.add_chromosome(1, endPos, chrName)

            if recType == RecordType.PSEUDOGENE:
                gene_obj: gm.Gene = gm.PseudoGene(
                    gid, note, startPos, endPos, chrName, strand, name, annots
                )
            else:
                gene_obj = gm.Gene(gid, note, startPos, endPos, chrName, strand, name, annots)

            try:
                other = model.allGenes[str(gene_obj)]
                conditionalException(
                    f"line {lineCtr}: gene {gene_obj.id} associated with multiple loci: "
                    f"{other.minpos}-{other.maxpos} and {startPos}-{endPos}"
                )
            except KeyError:
                pass

            model.allChr[chrName].update(gene_obj)
            model.add_gene(gene_obj)
            geneAlias[gene_obj.name.upper()] = gene_obj.id.upper()
            geneAlias[gene_obj.id.upper()] = gene_obj.name.upper()
            geneCount += 1

        elif recType in [RecordType.EXON, RecordType.PSEUDOGENIC_EXON]:
            if chrName not in model.model:
                continue

            parent_record: gm.Gene | gm.mRNA | None = None
            tried = set()
            for key in gm.POSSIBLE_GENE_FIELDS:
                if (key not in annots) or (annots[key] in tried):
                    continue
                isoName = annots[key]
                parent_record = model.get_parent(isoName, chrName)
                if parent_record:
                    break
                tried.add(isoName)
            if not parent_record:
                continue

            if parent_record.featureType == RecordType.MRNA:
                gene_obj = parent_record.parent
                if gene_obj is None:
                    conditionalException(
                        f"line {lineCtr}: mRNA parent is missing gene for exon record"
                    )
                    continue
            else:
                if not isinstance(parent_record, gm.Gene):
                    conditionalException(
                        f"line {lineCtr}: exon parent {parent_record} is not a gene"
                    )
                    continue
                gene_obj = parent_record

            isoform = None
            isoName = ""
            tried = set()
            for key in gm.POSSIBLE_FORM_FIELDS:
                if key not in annots:
                    continue
                if annots[key] in tried:
                    continue
                isoName = annots[key]
                if isoName in gene_obj.isoforms:
                    isoform = gene_obj.isoforms[isoName]
                    break
                tried.add(isoName)

            if not (isoform or isoName):
                continue

            if strand in gm.VALID_STRANDS and strand != gene_obj.strand:
                conditionalException(
                    f"line {lineCtr}: exon strand ({strand}) != gene "
                    f"strand ({gene_obj.strand}) for {gene_obj.id}"
                )
            else:
                strand = gene_obj.strand

            if not isoform:
                isoAttr = {
                    gm.PARENT_FIELD: gene_obj.id,
                    gm.NAME_FIELD: isoName,
                    gm.ID_FIELD: isoName,
                }
                isoform = gm.Isoform(isoName, startPos, endPos, chrName, strand, attr=isoAttr)
                isoCount += 1

            exon = gm.Exon(startPos, endPos, chrName, strand, annots)
            exonCount = exonCount + 1 if gene_obj.addExon(isoform, exon) else exonCount

        elif recType in [RecordType.MRNA, RecordType.PSEUDOGENIC_TRANSCRIPT]:
            if chrName not in model.model:
                conditionalException(
                    f"line {lineCtr}: mRNA with missing chromosome dictionary {chrName} "
                    f"(known: {','.join(model.model.keys())})"
                )

            if gm.ID_FIELD not in annots:
                conditionalException(f"line {lineCtr}: mRNA with missing ID")

            transcript_id = annots[gm.ID_FIELD].upper()
            parent_id = model.get_annotation(gm.PARENT_FIELD, annots)
            if parent_id is None:
                continue

            parent_id_upper = parent_id.upper()
            mrna_gene: gm.Gene | None = None
            parent_candidate = model.get_parent(parent_id_upper, chrName)
            if isinstance(parent_candidate, gm.Gene):
                mrna_gene = parent_candidate
            if not mrna_gene:
                alias = parent_id_upper
                try:
                    alias = geneAlias[parent_id_upper]
                    alias_candidate = model.get_parent(alias, chrName)
                    if isinstance(alias_candidate, gm.Gene):
                        mrna_gene = alias_candidate
                except KeyError:
                    if verbose:
                        if alias == parent_id_upper:
                            sys.stderr.write(
                                f"line {lineCtr}: no gene {parent_id_upper} found for {recType}\n"
                            )
                        else:
                            sys.stderr.write(
                                f"line {lineCtr}: no gene '{parent_id_upper}' or "
                                f"'{alias}' found for {recType}\n"
                            )
                    continue

            if mrna_gene is None:
                conditionalException(
                    f"line {lineCtr}: no gene parent found for transcript {parent_id_upper}"
                )
                continue

            if strand in gm.VALID_STRANDS and strand != mrna_gene.strand:
                conditionalException(
                    f"line {lineCtr}: mRNA strand ({strand}) does not "
                    f"match gene strand ({mrna_gene.strand})"
                )
            else:
                strand = mrna_gene.strand

            mrnaAttr = {
                gm.PARENT_FIELD: mrna_gene.id,
                gm.NAME_FIELD: transcript_id,
                gm.ID_FIELD: transcript_id,
            }
            mrna = gm.mRNA(transcript_id, startPos, endPos, chrName, strand, attr=mrnaAttr)
            mrna_gene.addmRNA(mrna)
            model.mRNAforms[chrName] = model.mRNAforms.setdefault(chrName, {})
            model.mRNAforms[chrName][transcript_id] = mrna
            mrnaCount += 1

        elif recType in gm.CDS_TYPES:
            if chrName not in model.model:
                conditionalException(
                    f"line {lineCtr}: {recType} has unrecognized chromosome: {chrName} "
                    f"(known: {','.join(model.model.keys())})"
                )
            if chrName not in model.mRNAforms:
                conditionalException(
                    f"line {lineCtr}: {recType} has unrecognized chromosome: {chrName} "
                    f"(known: {','.join(model.mRNAforms.keys())})"
                )

            mrna_record = model.get_parent(annots[gm.PARENT_FIELD], chrName, searchGenes=False)
            if not mrna_record:
                if verbose:
                    sys.stderr.write(
                        f"line {lineCtr}: no mRNA {annots[gm.PARENT_FIELD]} found for {recType}\n"
                    )
                continue
            if mrna_record.featureType != RecordType.MRNA:
                conditionalException(
                    f"line {lineCtr}: parent {annots[gm.PARENT_FIELD]} is not an mRNA record"
                )
                continue
            gene_obj = mrna_record.parent
            if gene_obj is None:
                conditionalException(
                    f"line {lineCtr}: mRNA {mrna_record.id} is missing a parent gene for CDS record"
                )
                continue

            if strand in gm.VALID_STRANDS and strand != mrna_record.strand:
                conditionalException(
                    f"line {lineCtr}: CDS strand ({strand}) does not "
                    f"match mRNA strand ({mrna_record.strand})"
                )
            else:
                strand = mrna_record.strand

            cds = gm.cdsFactory(recType, startPos, endPos, chrName, strand, annots)
            if gene_obj.addCDS(mrna_record, cds):
                cdsCount += 1

        elif recType == RecordType.CHROMOSOME:
            model.add_chromosome(startPos, endPos, chrName)

        elif recType in [RecordType.FIVE_PRIME_UTR, RecordType.THREE_PRIME_UTR]:
            if chrName not in model.model:
                continue
            if gm.PARENT_FIELD not in annots:
                continue

            parent_id = annots[gm.PARENT_FIELD]
            parent_record = model.get_parent(parent_id, chrName)
            if parent_record:
                if strand in gm.VALID_STRANDS and strand != parent_record.strand:
                    conditionalException(
                        f"line {lineCtr}: {recType} strand ({strand}) does not match "
                        f"parent strand ({parent_record.strand})"
                    )
                else:
                    strand = parent_record.strand
                parent_record.addFeature(
                    gm.BaseFeature(recType, startPos, endPos, chrName, strand, annots)
                )
            else:
                if verbose:
                    sys.stderr.write(f"line {lineCtr}: no parent {parent_id} found for {recType}\n")
                continue

        elif recType not in gm.IGNORE_RECTYPES:
            if chrName not in model.model:
                continue
            if gm.PARENT_FIELD not in annots:
                continue

            parent_gene_id = annots[gm.PARENT_FIELD].split(".")[0]
            if parent_gene_id not in model.model[chrName]:
                continue

            feature_gene = model.model[chrName][parent_gene_id]
            try:
                feature_gene.addFeature(
                    gm.BaseFeature(recType, startPos, endPos, chrName, strand, annots)
                )
            except ValueError as err:
                conditionalException(f"line {lineCtr}: {err}")

    indicator.finish()
    if verbose:
        if geneCount > 0:
            sys.stderr.write(
                "Loaded {} genes with {} isoforms, {} exons (avg. {:.1f}/gene), "
                "{} mRNA, {} CDS (avg. {:.1f}/gene)\n".format(
                    comma_format(geneCount),
                    comma_format(isoCount),
                    comma_format(exonCount),
                    float(exonCount) / geneCount,
                    comma_format(mrnaCount),
                    comma_format(cdsCount),
                    float(cdsCount / geneCount),
                )
            )
        else:
            sys.stderr.write("** Warning: no genes loaded!\n")


__all__ = ["load_gene_model_records"]
