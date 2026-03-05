"""Benchmark helpers for optional Polars GFF ingestion workflows."""

from __future__ import annotations

import hashlib
import json
import tracemalloc
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import TypedDict

from SpliceGrapher.formats.gene_model import GeneModel
from SpliceGrapher.formats.polars_gff import (
    PolarsNotInstalledError,
    load_gff_rows,
    load_gff_to_polars,
)

DEFAULT_DATASET_SIZES: dict[str, int] = {
    "small": 100,
    "medium": 750,
    "large": 3000,
}

RUNTIME_SPEEDUP_THRESHOLD = 0.20
MEMORY_REGRESSION_THRESHOLD = 0.10
REQUIRED_REAL_DATASETS = 3
MIN_WINNING_DATASETS = 2


@dataclass(frozen=True, slots=True)
class BenchmarkMetrics:
    """Single loader benchmark metrics."""

    mean_seconds: float
    max_seconds: float
    peak_mebibytes: float
    rows: int
    status: str = "ok"
    analytics_signature: str | None = None


@dataclass(frozen=True, slots=True)
class GoNoGoDecision:
    """Explicit decision record for the Polars go/no-go policy."""

    recommendation: str
    decision: str
    evaluated_real_datasets: int
    compared_real_datasets: int
    winning_datasets: int
    runtime_speedup_threshold: float = RUNTIME_SPEEDUP_THRESHOLD
    memory_regression_threshold: float = MEMORY_REGRESSION_THRESHOLD
    required_real_datasets: int = REQUIRED_REAL_DATASETS
    min_winning_datasets: int = MIN_WINNING_DATASETS
    notes: str = ""


@dataclass(frozen=True, slots=True)
class SingleCycleEvaluation:
    """Single-cycle benchmark output for synthetic + real workloads."""

    synthetic: dict[str, dict[str, dict[str, BenchmarkMetrics]]]
    real: dict[str, dict[str, dict[str, BenchmarkMetrics]]]
    decision: GoNoGoDecision


class TranscriptExonPayload(TypedDict):
    transcript: str
    strand: str
    chrom: str
    exons: list[tuple[int, int]]


class AnalyticsPayload(TypedDict):
    shared_exons: list[tuple[str, str, int, int, int]]
    transcript_exons: list[TranscriptExonPayload]
    junctions: list[tuple[str, int, int]]


def _as_int(value: int | str | None, *, field: str) -> int:
    if value is None:
        raise ValueError(f"Missing numeric value for {field}")
    return int(value)


def _count_data_rows(path: Path) -> int:
    with path.open("r", encoding="utf-8") as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith("#"))


def _measure_loader(
    loader: Callable[[Path], object],
    path: Path,
    *,
    iterations: int,
    rows: int,
) -> BenchmarkMetrics:
    # Warm up once so one-time import/init overhead does not dominate metrics.
    _ = loader(path)

    timings: list[float] = []
    peaks: list[int] = []

    for _ in range(iterations):
        tracemalloc.start()
        start = perf_counter()
        _ = loader(path)
        elapsed = perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        timings.append(elapsed)
        peaks.append(peak)

    return BenchmarkMetrics(
        mean_seconds=sum(timings) / len(timings),
        max_seconds=max(timings),
        peak_mebibytes=max(peaks) / (1024.0 * 1024.0),
        rows=rows,
    )


def _error_metrics(
    rows: int,
    *,
    exception: Exception,
    analytics_signature: str | None = None,
) -> BenchmarkMetrics:
    return BenchmarkMetrics(
        mean_seconds=0.0,
        max_seconds=0.0,
        peak_mebibytes=0.0,
        rows=rows,
        status=f"error:{type(exception).__name__}",
        analytics_signature=analytics_signature,
    )


def _measure_analytics_workload(
    workload: Callable[[Path], tuple[int, str]],
    path: Path,
    *,
    iterations: int,
) -> BenchmarkMetrics:
    # Warm up once so one-time import/init overhead does not dominate metrics.
    warm_rows, warm_signature = workload(path)

    timings: list[float] = []
    peaks: list[int] = []
    rows = warm_rows
    signature = warm_signature

    for _ in range(iterations):
        tracemalloc.start()
        start = perf_counter()
        rows_i, signature_i = workload(path)
        elapsed = perf_counter() - start
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()

        if signature_i != signature:
            raise ValueError("analytics signature drift within a single workload")

        rows = rows_i
        timings.append(elapsed)
        peaks.append(peak)

    return BenchmarkMetrics(
        mean_seconds=sum(timings) / len(timings),
        max_seconds=max(timings),
        peak_mebibytes=max(peaks) / (1024.0 * 1024.0),
        rows=rows,
        analytics_signature=signature,
    )


def _normalize_parent_ids(value: str | None) -> list[str]:
    if not value:
        return []
    return [part for part in value.split(",") if part]


def _exon_records_from_rows(path: Path) -> list[tuple[str, str, int, int, str]]:
    rows = load_gff_rows(path, ignore_malformed=False)
    result: list[tuple[str, str, int, int, str]] = []
    for row in rows:
        if row["type"] != "exon":
            continue
        chrom = str(row["chrom"])
        strand = str(row["strand"])
        start = _as_int(row["start"], field="start")
        end = _as_int(row["end"], field="end")
        parent_value = row["parent_id"] if isinstance(row["parent_id"], str) else None
        for parent_id in _normalize_parent_ids(parent_value):
            result.append((chrom, strand, start, end, parent_id))
    return result


def _exon_records_from_polars(path: Path) -> list[tuple[str, str, int, int, str]]:
    data_frame = load_gff_to_polars(path, ignore_malformed=False)
    result: list[tuple[str, str, int, int, str]] = []
    for row in data_frame.to_dicts():
        if row["type"] != "exon":
            continue
        chrom = str(row["chrom"])
        strand = str(row["strand"])
        start = int(row["start"])
        end = int(row["end"])
        for parent_id in _normalize_parent_ids(
            row["parent_id"] if isinstance(row["parent_id"], str) else None
        ):
            result.append((chrom, strand, start, end, parent_id))
    return result


def _exon_records_from_gene_model(path: Path) -> list[tuple[str, str, int, int, str]]:
    model = GeneModel(str(path), verbose=False)
    result: list[tuple[str, str, int, int, str]] = []

    for gene in model.all_genes.values():
        for isoform in gene.isoforms.values():
            for exon in isoform.exons:
                result.append((exon.chromosome, exon.strand, exon.minpos, exon.maxpos, isoform.id))

    return result


def _analytics_signature(exon_records: list[tuple[str, str, int, int, str]]) -> str:
    if not exon_records:
        empty_payload: AnalyticsPayload = {
            "shared_exons": [],
            "transcript_exons": [],
            "junctions": [],
        }
        encoded = json.dumps(empty_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    exon_group_map: dict[tuple[str, str, int, int], set[str]] = {}
    transcript_map: dict[str, list[tuple[int, int, str, str]]] = {}

    for chrom, strand, start, end, parent_id in exon_records:
        exon_group_map.setdefault((chrom, strand, start, end), set()).add(parent_id)
        transcript_map.setdefault(parent_id, []).append((start, end, strand, chrom))

    shared_exons = sorted(
        (chrom, strand, start, end, len(parents))
        for (chrom, strand, start, end), parents in exon_group_map.items()
        if len(parents) > 1
    )

    transcript_exons: list[TranscriptExonPayload] = []
    junctions: list[tuple[str, int, int]] = []

    for transcript_id, exons in sorted(transcript_map.items()):
        strand = exons[0][2]
        chrom = exons[0][3]
        ordered = sorted(((start, end) for start, end, _, _ in exons), reverse=(strand == "-"))
        transcript_exons.append(
            {
                "transcript": transcript_id,
                "strand": strand,
                "chrom": chrom,
                "exons": ordered,
            }
        )

        if len(ordered) < 2:
            continue

        prev_start, prev_end = ordered[0]
        for curr_start, curr_end in ordered[1:]:
            if strand == "+":
                donor = prev_end
                acceptor = curr_start - 2
            else:
                donor = prev_start - 2
                acceptor = curr_end
            junctions.append((transcript_id, donor, acceptor))
            prev_start, prev_end = curr_start, curr_end

    signature_payload: AnalyticsPayload = {
        "shared_exons": shared_exons,
        "transcript_exons": transcript_exons,
        "junctions": sorted(junctions),
    }
    encoded = json.dumps(signature_payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def _workload_rows(path: Path) -> tuple[int, str]:
    records = _exon_records_from_rows(path)
    return len(records), _analytics_signature(records)


def _workload_polars(path: Path) -> tuple[int, str]:
    records = _exon_records_from_polars(path)
    return len(records), _analytics_signature(records)


def _workload_gene_model(path: Path) -> tuple[int, str]:
    records = _exon_records_from_gene_model(path)
    return len(records), _analytics_signature(records)


def _benchmark_classes_for_path(
    path: Path,
    *,
    iterations: int,
    include_polars: bool,
) -> dict[str, dict[str, BenchmarkMetrics]]:
    return {
        "ingest": benchmark_gff_path(path, iterations=iterations, include_polars=include_polars),
        "end_to_end": benchmark_end_to_end_gff_path(
            path,
            iterations=iterations,
            include_polars=include_polars,
        ),
    }


def write_synthetic_gff(path: str | Path, *, gene_count: int, exons_per_gene: int = 3) -> int:
    """Write a deterministic synthetic GFF file and return record count."""
    if gene_count <= 0:
        raise ValueError("gene_count must be >= 1")
    if exons_per_gene <= 0:
        raise ValueError("exons_per_gene must be >= 1")

    out_path = Path(path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    record_count = 0
    with out_path.open("w", encoding="utf-8") as handle:
        handle.write("##gff-version 3\n")
        for gene_idx in range(1, gene_count + 1):
            gene_start = (gene_idx - 1) * 1000 + 1
            gene_end = gene_start + (exons_per_gene * 100)
            gene_id = f"GENE{gene_idx}"
            transcript_id = f"TX{gene_idx}"

            handle.write(
                f"chr1\tbenchmark\tgene\t{gene_start}\t{gene_end}\t.\t+\t."
                f"\tID={gene_id};Name={gene_id}\n"
            )
            record_count += 1

            handle.write(
                f"chr1\tbenchmark\tmrna\t{gene_start}\t{gene_end}\t.\t+\t."
                f"\tID={transcript_id};Parent={gene_id};Name={transcript_id}\n"
            )
            record_count += 1

            for exon_idx in range(1, exons_per_gene + 1):
                exon_start = gene_start + ((exon_idx - 1) * 100)
                exon_end = exon_start + 79
                exon_id = f"EX{gene_idx}.{exon_idx}"
                handle.write(
                    f"chr1\tbenchmark\texon\t{exon_start}\t{exon_end}\t.\t+\t."
                    f"\tID={exon_id};Parent={transcript_id}\n"
                )
                record_count += 1

    return record_count


def benchmark_gff_path(
    path: str | Path,
    *,
    iterations: int = 3,
    include_polars: bool = True,
) -> dict[str, BenchmarkMetrics]:
    """Benchmark available ingest loaders for a single GFF file."""
    if iterations <= 0:
        raise ValueError("iterations must be >= 1")

    gff_path = Path(path)
    rows = _count_data_rows(gff_path)

    results: dict[str, BenchmarkMetrics] = {}

    try:
        results["gene_model"] = _measure_loader(
            lambda p: GeneModel(str(p), verbose=False),
            gff_path,
            iterations=iterations,
            rows=rows,
        )
    except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
        results["gene_model"] = _error_metrics(rows, exception=exc)

    try:
        results["rows"] = _measure_loader(
            lambda p: load_gff_rows(p, ignore_malformed=False),
            gff_path,
            iterations=iterations,
            rows=rows,
        )
    except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
        results["rows"] = _error_metrics(rows, exception=exc)

    if include_polars:
        try:
            results["polars_df"] = _measure_loader(
                lambda p: load_gff_to_polars(p, ignore_malformed=False),
                gff_path,
                iterations=iterations,
                rows=rows,
            )
        except PolarsNotInstalledError:
            results["polars_df"] = BenchmarkMetrics(
                mean_seconds=0.0,
                max_seconds=0.0,
                peak_mebibytes=0.0,
                rows=rows,
                status="unavailable",
            )
        except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
            results["polars_df"] = _error_metrics(rows, exception=exc)

    return results


def benchmark_end_to_end_gff_path(
    path: str | Path,
    *,
    iterations: int = 3,
    include_polars: bool = True,
) -> dict[str, BenchmarkMetrics]:
    """Benchmark end-to-end analytics workloads for a single GFF file."""
    if iterations <= 0:
        raise ValueError("iterations must be >= 1")

    gff_path = Path(path)

    results: dict[str, BenchmarkMetrics] = {}

    try:
        results["gene_model"] = _measure_analytics_workload(
            _workload_gene_model,
            gff_path,
            iterations=iterations,
        )
    except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
        results["gene_model"] = _error_metrics(_count_data_rows(gff_path), exception=exc)

    try:
        results["rows"] = _measure_analytics_workload(
            _workload_rows,
            gff_path,
            iterations=iterations,
        )
    except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
        results["rows"] = _error_metrics(_count_data_rows(gff_path), exception=exc)

    if include_polars:
        try:
            results["polars_df"] = _measure_analytics_workload(
                _workload_polars,
                gff_path,
                iterations=iterations,
            )
        except PolarsNotInstalledError:
            results["polars_df"] = BenchmarkMetrics(
                mean_seconds=0.0,
                max_seconds=0.0,
                peak_mebibytes=0.0,
                rows=results["rows"].rows,
                status="unavailable",
                analytics_signature=None,
            )
        except Exception as exc:  # pragma: no cover - defensive for heterogeneous fixtures
            results["polars_df"] = _error_metrics(
                results["rows"].rows,
                exception=exc,
            )

    return results


def benchmark_matrix(
    work_dir: str | Path,
    *,
    dataset_sizes: dict[str, int] | None = None,
    exons_per_gene: int = 3,
    iterations: int = 3,
    include_polars: bool = True,
) -> dict[str, dict[str, BenchmarkMetrics]]:
    """Benchmark ingest loaders for multiple synthetic dataset sizes."""
    base_dir = Path(work_dir)
    base_dir.mkdir(parents=True, exist_ok=True)

    sizes = dataset_sizes if dataset_sizes is not None else DEFAULT_DATASET_SIZES
    result: dict[str, dict[str, BenchmarkMetrics]] = {}

    for label, gene_count in sizes.items():
        gff_path = base_dir / f"{label}.synthetic.gff3"
        write_synthetic_gff(gff_path, gene_count=gene_count, exons_per_gene=exons_per_gene)
        result[label] = benchmark_gff_path(
            gff_path,
            iterations=iterations,
            include_polars=include_polars,
        )

    return result


def evaluate_go_no_go(real_end_to_end: dict[str, dict[str, BenchmarkMetrics]]) -> GoNoGoDecision:
    """Apply explicit go/no-go rule for Polars adoption."""
    evaluated = len(real_end_to_end)
    if evaluated < REQUIRED_REAL_DATASETS:
        return GoNoGoDecision(
            recommendation="defer",
            decision="insufficient_data",
            evaluated_real_datasets=evaluated,
            compared_real_datasets=0,
            winning_datasets=0,
            notes="need at least three real datasets for final decision",
        )

    compared = 0
    winners = 0

    for loaders in real_end_to_end.values():
        baseline = loaders.get("rows")
        polars = loaders.get("polars_df")

        if not baseline or not polars:
            continue
        if baseline.status != "ok" or polars.status != "ok":
            continue

        compared += 1

        if baseline.analytics_signature != polars.analytics_signature:
            return GoNoGoDecision(
                recommendation="defer",
                decision="correctness_drift",
                evaluated_real_datasets=evaluated,
                compared_real_datasets=compared,
                winning_datasets=winners,
                notes="analytics signature mismatch between baseline and polars",
            )

        if baseline.mean_seconds <= 0:
            continue

        runtime_speedup = 1.0 - (polars.mean_seconds / baseline.mean_seconds)
        memory_ratio = (
            polars.peak_mebibytes / baseline.peak_mebibytes if baseline.peak_mebibytes > 0 else 1.0
        )

        if runtime_speedup >= RUNTIME_SPEEDUP_THRESHOLD and memory_ratio <= (
            1.0 + MEMORY_REGRESSION_THRESHOLD
        ):
            winners += 1

    if compared < REQUIRED_REAL_DATASETS:
        return GoNoGoDecision(
            recommendation="defer",
            decision="insufficient_data",
            evaluated_real_datasets=evaluated,
            compared_real_datasets=compared,
            winning_datasets=winners,
            notes="fewer than three comparable real datasets",
        )

    if winners >= MIN_WINNING_DATASETS:
        return GoNoGoDecision(
            recommendation="adopt",
            decision="threshold_met",
            evaluated_real_datasets=evaluated,
            compared_real_datasets=compared,
            winning_datasets=winners,
            notes="runtime and memory thresholds met on required dataset count",
        )

    return GoNoGoDecision(
        recommendation="defer",
        decision="threshold_not_met",
        evaluated_real_datasets=evaluated,
        compared_real_datasets=compared,
        winning_datasets=winners,
        notes="required runtime and memory win thresholds not met",
    )


def run_single_cycle_evaluation(
    *,
    synthetic_work_dir: str | Path,
    real_datasets: dict[str, str | Path],
    synthetic_dataset_sizes: dict[str, int] | None = None,
    exons_per_gene: int = 3,
    iterations: int = 3,
    include_polars: bool = True,
) -> SingleCycleEvaluation:
    """Run one complete synthetic + real benchmark cycle."""
    work_dir = Path(synthetic_work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)

    sizes = (
        synthetic_dataset_sizes if synthetic_dataset_sizes is not None else DEFAULT_DATASET_SIZES
    )

    synthetic: dict[str, dict[str, dict[str, BenchmarkMetrics]]] = {}
    for label, gene_count in sizes.items():
        path = work_dir / f"{label}.synthetic.gff3"
        write_synthetic_gff(path, gene_count=gene_count, exons_per_gene=exons_per_gene)
        synthetic[label] = _benchmark_classes_for_path(
            path,
            iterations=iterations,
            include_polars=include_polars,
        )

    real: dict[str, dict[str, dict[str, BenchmarkMetrics]]] = {}
    for label, path_value in real_datasets.items():
        path = Path(path_value)
        if not path.is_file():
            raise FileNotFoundError(f"real dataset path not found: {path}")
        real[label] = _benchmark_classes_for_path(
            path,
            iterations=iterations,
            include_polars=include_polars,
        )

    end_to_end_real = {name: classes["end_to_end"] for name, classes in real.items()}
    decision = evaluate_go_no_go(end_to_end_real)

    return SingleCycleEvaluation(synthetic=synthetic, real=real, decision=decision)


def _metrics_to_dict(metrics: BenchmarkMetrics) -> dict[str, object]:
    return {
        "status": metrics.status,
        "rows": metrics.rows,
        "mean_seconds": metrics.mean_seconds,
        "max_seconds": metrics.max_seconds,
        "peak_mebibytes": metrics.peak_mebibytes,
        "analytics_signature": metrics.analytics_signature,
    }


def evaluation_to_json_dict(evaluation: SingleCycleEvaluation) -> dict[str, object]:
    """Convert single-cycle evaluation to JSON-safe dictionary payload."""

    def convert_dataset(
        section: dict[str, dict[str, dict[str, BenchmarkMetrics]]],
    ) -> dict[str, dict[str, dict[str, dict[str, object]]]]:
        payload: dict[str, dict[str, dict[str, dict[str, object]]]] = {}
        for dataset, classes in section.items():
            payload[dataset] = {}
            for class_name, loaders in classes.items():
                payload[dataset][class_name] = {
                    loader_name: _metrics_to_dict(loader_metrics)
                    for loader_name, loader_metrics in loaders.items()
                }
        return payload

    return {
        "synthetic": convert_dataset(evaluation.synthetic),
        "real": convert_dataset(evaluation.real),
        "decision": {
            "recommendation": evaluation.decision.recommendation,
            "decision": evaluation.decision.decision,
            "evaluated_real_datasets": evaluation.decision.evaluated_real_datasets,
            "compared_real_datasets": evaluation.decision.compared_real_datasets,
            "winning_datasets": evaluation.decision.winning_datasets,
            "runtime_speedup_threshold": evaluation.decision.runtime_speedup_threshold,
            "memory_regression_threshold": evaluation.decision.memory_regression_threshold,
            "required_real_datasets": evaluation.decision.required_real_datasets,
            "min_winning_datasets": evaluation.decision.min_winning_datasets,
            "notes": evaluation.decision.notes,
        },
    }


def matrix_to_markdown(matrix: dict[str, dict[str, BenchmarkMetrics]]) -> str:
    """Render benchmark matrix as a markdown table."""
    lines = [
        "| dataset | loader | status | rows | mean_s | max_s | peak_mib |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]

    for dataset, loaders in matrix.items():
        for loader_name, metrics in loaders.items():
            lines.append(
                "| "
                f"{dataset} | {loader_name} | {metrics.status} | {metrics.rows} | "
                f"{metrics.mean_seconds:.6f} | {metrics.max_seconds:.6f} | "
                f"{metrics.peak_mebibytes:.3f} |"
            )

    return "\n".join(lines)


def _section_to_markdown(
    title: str,
    section: dict[str, dict[str, dict[str, BenchmarkMetrics]]],
) -> str:
    blocks = [f"## {title}"]
    for dataset, classes in section.items():
        blocks.append(f"### {dataset}")
        for class_name, loaders in classes.items():
            blocks.append(f"#### {class_name}")
            blocks.append(matrix_to_markdown({dataset: loaders}))
            blocks.append("")
    return "\n".join(blocks).rstrip()


def evaluation_to_markdown(evaluation: SingleCycleEvaluation) -> str:
    """Render single-cycle evaluation as markdown report."""
    decision = evaluation.decision
    lines = [
        "# Polars GFF Single-Cycle Evaluation",
        "",
        "## Decision",
        f"- recommendation: `{decision.recommendation}`",
        f"- decision: `{decision.decision}`",
        f"- evaluated real datasets: `{decision.evaluated_real_datasets}`",
        f"- compared real datasets: `{decision.compared_real_datasets}`",
        f"- winning datasets: `{decision.winning_datasets}`",
        f"- thresholds: runtime>={int(decision.runtime_speedup_threshold * 100)}%, "
        f"memory<=+{int(decision.memory_regression_threshold * 100)}%",
    ]
    if decision.notes:
        lines.append(f"- notes: {decision.notes}")

    lines.extend(
        [
            "",
            _section_to_markdown("Synthetic", evaluation.synthetic),
            "",
            _section_to_markdown("Real", evaluation.real),
        ]
    )
    return "\n".join(lines)


__all__ = [
    "DEFAULT_DATASET_SIZES",
    "MEMORY_REGRESSION_THRESHOLD",
    "MIN_WINNING_DATASETS",
    "REQUIRED_REAL_DATASETS",
    "RUNTIME_SPEEDUP_THRESHOLD",
    "BenchmarkMetrics",
    "GoNoGoDecision",
    "SingleCycleEvaluation",
    "benchmark_end_to_end_gff_path",
    "benchmark_gff_path",
    "benchmark_matrix",
    "evaluate_go_no_go",
    "evaluation_to_json_dict",
    "evaluation_to_markdown",
    "matrix_to_markdown",
    "run_single_cycle_evaluation",
    "write_synthetic_gff",
]
