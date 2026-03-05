from __future__ import annotations

import hashlib
import tempfile
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path
from typing import TypedDict

import gffutils

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import AttrKey, RecordType
from SpliceGrapher.formats.gene_model import (
    CDS,
    FP_UTR,
    TP_UTR,
    Exon,
    Gene,
    GeneModel,
    Isoform,
    featureSortKey,
    mRNA,
)

GENE_FEATURE_TYPES = {RecordType.GENE, RecordType.PREDICTED_GENE, RecordType.PSEUDOGENE}
TRANSCRIPT_FEATURE_TYPES = {RecordType.MRNA, RecordType.PSEUDOGENIC_TRANSCRIPT}
CDS_FEATURE_TYPES = {
    RecordType.CDS,
    RecordType.FIVE_PRIME_UTR,
    RecordType.THREE_PRIME_UTR,
}

# External annotation schema keys (GFF/GTF), centralized to avoid literal drift.
TRANSCRIPT_ID_KEY = "transcript_id"
GENE_ID_KEY = "gene_id"
GENE_KEY = "gene"
GENE_NAME_KEY = "gene_name"
GENE_BIOTYPE_KEY = "gene_biotype"

TRANSCRIPT_ALIAS = "transcript"

TRANSCRIPT_ID_KEYS = (TRANSCRIPT_ID_KEY, AttrKey.ID, AttrKey.NAME)
GENE_ID_KEYS = (GENE_ID_KEY, AttrKey.ID, AttrKey.NAME, GENE_KEY, GENE_NAME_KEY)
GENE_NAME_KEYS = (AttrKey.NAME, GENE_NAME_KEY, GENE_ID_KEY)
GENE_NOTE_KEYS = (AttrKey.NOTE, GENE_BIOTYPE_KEY)
GENE_REFERENCE_KEYS = (GENE_ID_KEY, GENE_KEY, GENE_NAME_KEY)
TRANSCRIPT_OR_PARENT_KEYS = (TRANSCRIPT_ID_KEY, AttrKey.PARENT, AttrKey.ID)


class GeneRecord(TypedDict):
    id: str
    chrom: str
    strand: str
    name: str
    note: str | None
    minpos: int
    maxpos: int
    attrs: dict[str, str]


def _normalize_feature_type(raw_feature_type: str) -> RecordType | None:
    """Normalize raw feature type strings to canonical ``RecordType`` values."""
    normalized = raw_feature_type.strip().casefold()
    if not normalized:
        return None
    if normalized == TRANSCRIPT_ALIAS:
        return RecordType.MRNA
    try:
        return coerce_enum(normalized, RecordType, field="feature_type")
    except ValueError:
        return None


def _first_attr(feature: gffutils.Feature, keys: Sequence[str]) -> str | None:
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
        feature_type = _normalize_feature_type(feature.featuretype)
        if feature_type not in TRANSCRIPT_FEATURE_TYPES:
            continue

        transcript_id = _first_attr(feature, TRANSCRIPT_ID_KEYS) or feature.id
        if not transcript_id:
            continue

        parents = feature.attributes.get(AttrKey.PARENT, [])
        gene_id: str | None = parents[0] if parents else _first_attr(feature, GENE_REFERENCE_KEYS)
        if not gene_id:
            continue

        mapping[transcript_id] = gene_id
    return mapping


def _extract_gene_records(db: gffutils.FeatureDB) -> dict[str, GeneRecord]:
    """Collect normalized gene metadata from gene-like annotation records."""
    records: dict[str, GeneRecord] = {}
    for feature in db.all_features(order_by=("seqid", "start", "end")):
        feature_type = _normalize_feature_type(feature.featuretype)
        if feature_type not in GENE_FEATURE_TYPES:
            continue

        raw_gene_id = _first_attr(feature, GENE_ID_KEYS) or feature.id
        if raw_gene_id is None:
            continue

        gene_id = raw_gene_id.upper()
        chrom = feature.seqid.lower()
        strand = feature.strand if feature.strand in {"+", "-", "."} else "."
        attrs = _feature_attr_map(feature)
        name = _first_attr(feature, GENE_NAME_KEYS) or raw_gene_id
        note = _first_attr(feature, GENE_NOTE_KEYS)

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
    gene_records: dict[str, GeneRecord],
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

    rec: GeneRecord | None = gene_records.get(gene_id)
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
        rec["minpos"],
        rec["maxpos"],
        rec["chrom"],
        rec["strand"],
        rec["name"],
        rec["attrs"],
    )
    model.add_gene(gene)
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
        feature_type = _normalize_feature_type(feature.featuretype)

        if feature_type in TRANSCRIPT_FEATURE_TYPES:
            transcript_id = _first_attr(feature, TRANSCRIPT_ID_KEYS) or feature.id
            if transcript_id:
                transcript_meta[transcript_id] = (chrom, strand)
                maybe_gene = _first_attr(feature, GENE_REFERENCE_KEYS)
                parents = feature.attributes.get(AttrKey.PARENT)
                if parents:
                    transcript_gene.setdefault(transcript_id, parents[0])
                elif maybe_gene:
                    transcript_gene.setdefault(transcript_id, maybe_gene)

        if feature_type == RecordType.EXON:
            transcript_id = _first_attr(feature, TRANSCRIPT_OR_PARENT_KEYS)
            if transcript_id is None:
                continue
            exon_groups[transcript_id].append(feature)
            transcript_meta.setdefault(transcript_id, (chrom, strand))
            maybe_gene = _first_attr(feature, GENE_REFERENCE_KEYS)
            parents = feature.attributes.get(AttrKey.PARENT)
            if maybe_gene:
                transcript_gene.setdefault(transcript_id, maybe_gene)
            elif parents and transcript_id not in transcript_gene:
                transcript_gene[transcript_id] = parents[0]

        if feature_type in CDS_FEATURE_TYPES:
            transcript_id = _first_attr(feature, TRANSCRIPT_OR_PARENT_KEYS)
            if transcript_id is None:
                continue
            cds_groups[transcript_id].append(feature)
            transcript_meta.setdefault(transcript_id, (chrom, strand))
            maybe_gene = _first_attr(feature, GENE_REFERENCE_KEYS)
            if maybe_gene:
                transcript_gene.setdefault(transcript_id, maybe_gene)

    for chrom, maxpos in chrom_max.items():
        model.add_chromosome(1, maxpos, chrom)

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
            raw_gene_id = _first_attr(exons[0], GENE_REFERENCE_KEYS) if exons else None
        if raw_gene_id is None:
            raw_gene_id = f"gene_{transcript_id}"
        gene_id = raw_gene_id.upper()

        gene = _get_or_create_gene(model, gene_records, gene_id, chrom, strand, minpos, maxpos)

        iso_attr = {
            AttrKey.PARENT: gene.id,
            AttrKey.NAME: transcript_id,
            AttrKey.ID: transcript_id,
        }
        isoform = Isoform(transcript_id, minpos, maxpos, chrom, strand, attr=iso_attr)
        for exon_feature in exons:
            exon = Exon(
                exon_feature.start,
                exon_feature.end,
                chrom,
                strand,
                _feature_attr_map(exon_feature),
            )
            gene.add_exon(isoform, exon)

        if cds_records:
            mrna_attr = {
                AttrKey.PARENT: gene.id,
                AttrKey.NAME: transcript_id,
                AttrKey.ID: transcript_id,
            }
            mrna_record = mRNA(transcript_id, minpos, maxpos, chrom, strand, attr=mrna_attr)
            for cds_feature in cds_records:
                feature_type = _normalize_feature_type(cds_feature.featuretype)
                if feature_type == RecordType.CDS:
                    cds = CDS(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                elif feature_type == RecordType.FIVE_PRIME_UTR:
                    cds = FP_UTR(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                elif feature_type == RecordType.THREE_PRIME_UTR:
                    cds = TP_UTR(
                        cds_feature.start,
                        cds_feature.end,
                        chrom,
                        strand,
                        _feature_attr_map(cds_feature),
                    )
                else:
                    continue
                gene.add_cds(mrna_record, cds)

    model.make_sorted_model()
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
            model.get_all_genes(),
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
