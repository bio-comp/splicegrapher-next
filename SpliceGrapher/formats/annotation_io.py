from __future__ import annotations

import hashlib
import tempfile
from collections import defaultdict
from pathlib import Path

import gffutils

from SpliceGrapher.formats.GeneModel import (
    CDS,
    EXON_TYPE,
    FP_UTR,
    GENE_TYPE,
    ID_FIELD,
    MRNA_TYPE,
    NAME_FIELD,
    PARENT_FIELD,
    TP_UTR,
    Exon,
    Gene,
    GeneModel,
    Isoform,
    featureSortKey,
    mRNA,
)

GENE_FEATURE_TYPES = {GENE_TYPE, "predicted_gene", "pseudogene"}
TRANSCRIPT_FEATURE_TYPES = {MRNA_TYPE, "transcript", "pseudogenic_transcript"}
CDS_FEATURE_TYPES = {"cds", "five_prime_utr", "three_prime_utr"}


def _first_attr(feature: gffutils.Feature, keys: list[str]) -> str | None:
    """Return the first populated attribute value from ``keys``."""
    for key in keys:
        values = feature.attributes.get(key)
        if values:
            return values[0]
    return None


def _feature_attr_map(feature: gffutils.Feature) -> dict[str, str]:
    """Flatten gffutils multi-value attributes to first-value strings."""
    return {key: values[0] for key, values in feature.attributes.items() if values}


def _stable_db_path(source_path: Path, cache_dir: Path | None) -> Path:
    """Compute a deterministic cache path for the gffutils sqlite database."""
    if cache_dir is None:
        tmpdir = Path(tempfile.mkdtemp(prefix="idiffir_gffutils_"))
        return tmpdir / "annotations.db"

    cache_dir.mkdir(parents=True, exist_ok=True)
    digest = hashlib.sha256(
        f"{source_path.resolve()}::{source_path.stat().st_mtime_ns}::{source_path.stat().st_size}".encode(
            "utf-8"
        )
    ).hexdigest()[:16]
    return cache_dir / f"annotation_{digest}.db"


def _create_db(path: Path, cache_dir: Path | None = None) -> gffutils.FeatureDB:
    """Build and return a fresh ``gffutils.FeatureDB`` from ``path``."""
    db_path = _stable_db_path(path, cache_dir)
    gffutils.create_db(
        str(path),
        dbfn=str(db_path),
        force=True,
        keep_order=True,
        merge_strategy="create_unique",
        sort_attribute_values=True,
        disable_infer_genes=False,
        disable_infer_transcripts=False,
    )
    return gffutils.FeatureDB(str(db_path), keep_order=True)


def _transcript_gene_map(db: gffutils.FeatureDB) -> dict[str, str]:
    """Map transcript IDs to gene IDs from transcript-like annotation records."""
    mapping: dict[str, str] = {}
    for feature in db.all_features(order_by=("seqid", "start", "end")):
        ftype = feature.featuretype.lower()
        if ftype not in TRANSCRIPT_FEATURE_TYPES:
            continue

        transcript_id = _first_attr(feature, ["transcript_id", "ID", "Name"]) or feature.id
        if not transcript_id:
            continue

        gene_id = None
        parents = feature.attributes.get("Parent")
        if parents:
            gene_id = parents[0]
        if gene_id is None:
            gene_id = _first_attr(feature, ["gene_id", "gene", "gene_name", "Name"])
        if gene_id is None:
            continue

        mapping[transcript_id] = gene_id
    return mapping


def _extract_gene_records(db: gffutils.FeatureDB) -> dict[str, dict[str, object]]:
    """Collect normalized gene metadata from gene-like annotation records."""
    records: dict[str, dict[str, object]] = {}
    for feature in db.all_features(order_by=("seqid", "start", "end")):
        ftype = feature.featuretype.lower()
        if ftype not in GENE_FEATURE_TYPES:
            continue

        raw_gene_id = (
            _first_attr(feature, ["gene_id", "ID", "Name", "gene", "gene_name"]) or feature.id
        )
        if raw_gene_id is None:
            continue

        gene_id = raw_gene_id.upper()
        chrom = feature.seqid.lower()
        strand = feature.strand if feature.strand in {"+", "-", "."} else "."
        attrs = _feature_attr_map(feature)
        name = _first_attr(feature, ["Name", "gene_name", "gene_id"]) or raw_gene_id
        note = _first_attr(feature, ["Note", "gene_biotype"])

        if gene_id not in records:
            records[gene_id] = {
                "id": gene_id,
                "chrom": chrom,
                "strand": strand,
                "name": name,
                "note": note,
                "minpos": feature.start,
                "maxpos": feature.end,
                "attrs": attrs,
            }
        else:
            records[gene_id]["minpos"] = min(int(records[gene_id]["minpos"]), feature.start)
            records[gene_id]["maxpos"] = max(int(records[gene_id]["maxpos"]), feature.end)
    return records


def _get_or_create_gene(
    model: GeneModel,
    gene_records: dict[str, dict[str, object]],
    gene_id: str,
    chrom: str,
    strand: str,
    default_start: int,
    default_end: int,
) -> Gene:
    """Return an existing gene entry or create it from normalized gene metadata."""
    if chrom not in model.model:
        model.model[chrom] = {}

    existing = model.model[chrom].get(gene_id)
    if existing is not None:
        existing.minpos = min(existing.minpos, default_start)
        existing.maxpos = max(existing.maxpos, default_end)
        return existing

    rec = gene_records.get(gene_id)
    if rec is None:
        rec = {
            "id": gene_id,
            "chrom": chrom,
            "strand": strand,
            "name": gene_id,
            "note": None,
            "minpos": default_start,
            "maxpos": default_end,
            "attrs": {},
        }

    gene = Gene(
        rec["id"],
        rec["note"],
        int(rec["minpos"]),
        int(rec["maxpos"]),
        rec["chrom"],
        rec["strand"],
        rec["name"],
        rec["attrs"],
    )
    model.addGene(gene)
    return gene


def _build_gene_model_from_db(db: gffutils.FeatureDB) -> GeneModel:
    """Convert a gffutils feature database into a legacy ``GeneModel`` object."""
    model = GeneModel(None)
    transcript_gene = _transcript_gene_map(db)
    gene_records = _extract_gene_records(db)

    exon_groups: dict[str, list[gffutils.Feature]] = defaultdict(list)
    cds_groups: dict[str, list[gffutils.Feature]] = defaultdict(list)
    chrom_max: dict[str, int] = {}
    transcript_meta: dict[str, tuple[str, str]] = {}

    for feature in db.all_features(order_by=("seqid", "start", "end")):
        chrom = feature.seqid.lower()
        chrom_max[chrom] = max(chrom_max.get(chrom, 0), feature.end)
        strand = feature.strand if feature.strand in {"+", "-", "."} else "."
        ftype = feature.featuretype.lower()

        if ftype in TRANSCRIPT_FEATURE_TYPES:
            transcript_id = _first_attr(feature, ["transcript_id", "ID", "Name"]) or feature.id
            if transcript_id:
                transcript_meta[transcript_id] = (chrom, strand)
                maybe_gene = _first_attr(feature, ["gene_id", "gene", "gene_name"])
                parents = feature.attributes.get("Parent")
                if parents:
                    transcript_gene.setdefault(transcript_id, parents[0])
                elif maybe_gene:
                    transcript_gene.setdefault(transcript_id, maybe_gene)

        if ftype == EXON_TYPE:
            transcript_id = _first_attr(feature, ["transcript_id", "Parent", "ID"])
            if transcript_id is None:
                continue
            exon_groups[transcript_id].append(feature)
            transcript_meta.setdefault(transcript_id, (chrom, strand))
            maybe_gene = _first_attr(feature, ["gene_id", "gene", "gene_name"])
            parents = feature.attributes.get("Parent")
            if maybe_gene:
                transcript_gene.setdefault(transcript_id, maybe_gene)
            elif parents and transcript_id not in transcript_gene:
                transcript_gene[transcript_id] = parents[0]

        if ftype in CDS_FEATURE_TYPES:
            transcript_id = _first_attr(feature, ["transcript_id", "Parent", "ID"])
            if transcript_id is None:
                continue
            cds_groups[transcript_id].append(feature)
            transcript_meta.setdefault(transcript_id, (chrom, strand))
            maybe_gene = _first_attr(feature, ["gene_id", "gene", "gene_name"])
            if maybe_gene:
                transcript_gene.setdefault(transcript_id, maybe_gene)

    for chrom, maxpos in chrom_max.items():
        model.addChromosome(1, maxpos, chrom)

    all_transcripts = sorted(
        set(exon_groups.keys()) | set(cds_groups.keys()) | set(transcript_meta.keys())
    )
    for transcript_id in all_transcripts:
        exons = sorted(exon_groups.get(transcript_id, []), key=lambda f: (f.start, f.end))
        cds_records = sorted(cds_groups.get(transcript_id, []), key=lambda f: (f.start, f.end))
        if not exons and not cds_records:
            continue

        chrom, strand = transcript_meta.get(transcript_id, (None, None))
        if chrom is None and exons:
            chrom = exons[0].seqid.lower()
        if strand is None and exons:
            strand = exons[0].strand if exons[0].strand in {"+", "-", "."} else "."
        if chrom is None:
            continue
        if strand is None:
            strand = "."

        starts = [f.start for f in exons] + [f.start for f in cds_records]
        ends = [f.end for f in exons] + [f.end for f in cds_records]
        minpos = min(starts)
        maxpos = max(ends)

        raw_gene_id = transcript_gene.get(transcript_id)
        if raw_gene_id is None:
            raw_gene_id = _first_attr(exons[0], ["gene_id", "gene", "gene_name"]) if exons else None
        if raw_gene_id is None:
            raw_gene_id = f"gene_{transcript_id}"
        gene_id = raw_gene_id.upper()

        gene = _get_or_create_gene(model, gene_records, gene_id, chrom, strand, minpos, maxpos)

        iso_attr = {PARENT_FIELD: gene.id, NAME_FIELD: transcript_id, ID_FIELD: transcript_id}
        isoform = Isoform(transcript_id, minpos, maxpos, chrom, strand, attr=iso_attr)
        for exon_feature in exons:
            exon = Exon(
                exon_feature.start,
                exon_feature.end,
                chrom,
                strand,
                _feature_attr_map(exon_feature),
            )
            gene.addExon(isoform, exon)

        if cds_records:
            mrna_attr = {PARENT_FIELD: gene.id, NAME_FIELD: transcript_id, ID_FIELD: transcript_id}
            mrna_record = mRNA(transcript_id, minpos, maxpos, chrom, strand, attr=mrna_attr)
            for cds_feature in cds_records:
                ftype = cds_feature.featuretype.lower()
                if ftype == "cds":
                    cds = CDS(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                elif ftype == "five_prime_utr":
                    cds = FP_UTR(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                else:
                    cds = TP_UTR(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                gene.addCDS(mrna_record, cds)

    model.makeSortedModel()
    return model


def _iter_transcript_exons(gene: Gene):
    """Yield transcript IDs and sorted exon lists for intron derivation."""
    for iso in gene.isoforms.values():
        exons = sorted(iso.exons, key=featureSortKey)
        if len(exons) > 1:
            yield iso.id, exons

    for mrna_rec in gene.mrna.values():
        exons = sorted(mrna_rec.sortedExons(), key=featureSortKey)
        if len(exons) > 1:
            yield mrna_rec.id, exons


def write_intron_cache(model: GeneModel, outdir: str | Path) -> Path:
    """Write a deterministic intron BED cache derived from model exon topology."""
    outdir_path = Path(outdir)
    cache_dir = outdir_path / ".idiffir_cache"
    cache_dir.mkdir(parents=True, exist_ok=True)
    cache_path = cache_dir / "idiffir_introns.bed"

    with cache_path.open("w", encoding="utf-8") as handle:
        for gene in sorted(
            model.getAllGenes(),
            key=lambda g: (g.chromosome, g.minpos, g.maxpos, g.id),
        ):
            for transcript_id, exons in _iter_transcript_exons(gene):
                for intron_index, (left, right) in enumerate(zip(exons[:-1], exons[1:]), start=1):
                    intron_start = left.maxpos + 1
                    intron_end = right.minpos - 1
                    if intron_start > intron_end:
                        continue
                    # BED uses 0-based, half-open coordinates.
                    handle.write(
                        f"{gene.chromosome}\t{intron_start - 1}\t{intron_end}\t"
                        f"{gene.id}\t{transcript_id}\t{intron_index}\t{gene.strand}\n"
                    )

    return cache_path


def load_gene_models(path: str, **args) -> GeneModel:
    """Load GTF/GFF annotations through gffutils and return a ``GeneModel``."""
    model_path = Path(path)
    if not model_path.exists():
        raise ValueError(f"Gene model file not found: {path}")

    cache_dir = args.get("cache_dir")
    db = _create_db(model_path, Path(cache_dir) if cache_dir else None)
    model = _build_gene_model_from_db(db)

    outdir = args.get("outdir")
    if outdir:
        write_intron_cache(model, outdir)

    return model
