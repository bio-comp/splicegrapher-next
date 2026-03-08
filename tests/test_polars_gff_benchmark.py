from __future__ import annotations

import json
import tomllib
from pathlib import Path
from typing import cast

from benchmarks import polars_gff_benchmark as benchmark
from SpliceGrapher.formats import polars_gff


def test_wheel_packages_only_runtime_namespace() -> None:
    pyproject_path = Path(__file__).resolve().parents[1] / "pyproject.toml"
    config = tomllib.loads(pyproject_path.read_text(encoding="utf-8"))

    wheel_packages = config["tool"]["hatch"]["build"]["targets"]["wheel"]["packages"]
    assert wheel_packages == ["SpliceGrapher"]


def test_write_synthetic_gff_generates_expected_record_count(tmp_path: Path) -> None:
    gff_path = tmp_path / "synthetic.gff3"
    row_count = benchmark.write_synthetic_gff(gff_path, gene_count=3, exons_per_gene=2)

    # per gene: gene + mrna + exons
    assert row_count == 3 * (1 + 1 + 2)

    rows = polars_gff.load_gff_rows(gff_path)
    assert len(rows) == row_count
    assert rows[0]["type"] == "gene"
    assert rows[1]["type"] == "mrna"


def test_benchmark_gff_path_without_polars_returns_core_loaders(tmp_path: Path) -> None:
    gff_path = tmp_path / "sample.gff3"
    benchmark.write_synthetic_gff(gff_path, gene_count=2, exons_per_gene=2)

    results = benchmark.benchmark_gff_path(gff_path, iterations=1, include_polars=False)

    assert set(results.keys()) == {"gene_model", "rows"}
    assert results["gene_model"].status == "ok"
    assert results["rows"].status == "ok"


def test_benchmark_matrix_marks_polars_unavailable(tmp_path: Path, monkeypatch) -> None:
    def _raise_missing(path: str | Path, *, ignore_malformed: bool = False):
        raise polars_gff.PolarsNotInstalledError("missing")

    monkeypatch.setattr(benchmark, "load_gff_to_polars", _raise_missing)

    matrix = benchmark.benchmark_matrix(
        tmp_path,
        dataset_sizes={"tiny": 2},
        exons_per_gene=2,
        iterations=1,
        include_polars=True,
    )

    tiny = matrix["tiny"]
    assert tiny["polars_df"].status == "unavailable"
    assert tiny["gene_model"].status == "ok"
    assert tiny["rows"].status == "ok"


def test_end_to_end_benchmark_exposes_analytics_signature(tmp_path: Path) -> None:
    gff_path = tmp_path / "sample.gff3"
    benchmark.write_synthetic_gff(gff_path, gene_count=3, exons_per_gene=2)

    results = benchmark.benchmark_end_to_end_gff_path(
        gff_path,
        iterations=1,
        include_polars=False,
    )

    assert set(results.keys()) == {"gene_model", "rows"}
    assert results["gene_model"].analytics_signature is not None
    assert results["rows"].analytics_signature is not None
    assert results["gene_model"].analytics_signature == results["rows"].analytics_signature


def test_go_no_go_requires_three_real_datasets() -> None:
    metrics = benchmark.BenchmarkMetrics(
        mean_seconds=0.1,
        max_seconds=0.1,
        peak_mebibytes=1.0,
        rows=100,
        status="ok",
        analytics_signature="abc123",
    )
    real_results = {
        "real1": {"rows": metrics, "polars_df": metrics},
        "real2": {"rows": metrics, "polars_df": metrics},
    }

    decision = benchmark.evaluate_go_no_go(real_results)
    assert decision.recommendation == "defer"
    assert decision.decision == "insufficient_data"


def test_evaluation_json_contains_decision_and_two_classes(tmp_path: Path) -> None:
    gff_path = tmp_path / "real.gff3"
    benchmark.write_synthetic_gff(gff_path, gene_count=2, exons_per_gene=2)

    evaluation = benchmark.run_single_cycle_evaluation(
        synthetic_work_dir=tmp_path / "synthetic",
        real_datasets={"realA": gff_path, "realB": gff_path, "realC": gff_path},
        synthetic_dataset_sizes={"small": 2},
        exons_per_gene=2,
        iterations=1,
        include_polars=False,
    )

    payload = benchmark.evaluation_to_json_dict(evaluation)
    synthetic_payload = cast(dict[str, dict[str, object]], payload["synthetic"])
    real_payload = cast(dict[str, dict[str, object]], payload["real"])
    decision_payload = cast(dict[str, object], payload["decision"])
    assert set(payload.keys()) == {"synthetic", "real", "decision"}
    assert set(synthetic_payload["small"].keys()) == {"ingest", "end_to_end"}
    assert set(real_payload["realA"].keys()) == {"ingest", "end_to_end"}
    assert isinstance(decision_payload["recommendation"], str)
    # Ensure payload is cleanly serializable for report emission.
    json.dumps(payload)
