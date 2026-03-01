from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

from SpliceGrapher.formats.polars_gff_benchmark import write_synthetic_gff


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _script_path() -> Path:
    return _repo_root() / "scripts" / "benchmarks" / "run_gff_loader_benchmarks.py"


def _create_real_dataset_files(tmp_path: Path, *, count: int) -> list[Path]:
    result: list[Path] = []
    for idx in range(1, count + 1):
        path = tmp_path / f"real{idx}.gff3"
        write_synthetic_gff(path, gene_count=idx + 1, exons_per_gene=2)
        result.append(path)
    return result


def test_runner_writes_outputs_with_skip_polars(tmp_path: Path) -> None:
    script = _script_path()
    synthetic_dir = tmp_path / "synthetic"
    json_out = tmp_path / "result.json"
    markdown_out = tmp_path / "result.md"
    real_paths = _create_real_dataset_files(tmp_path, count=3)

    cmd = [
        sys.executable,
        str(script),
        "--iterations",
        "1",
        "--synthetic-work-dir",
        str(synthetic_dir),
        "--size",
        "small=2",
        "--skip-polars",
        "--json-out",
        str(json_out),
        "--markdown-out",
        str(markdown_out),
    ]
    for idx, path in enumerate(real_paths, start=1):
        cmd.extend(["--real-dataset", f"real{idx}={path}"])

    result = subprocess.run(
        cmd,
        cwd=_repo_root(),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0
    assert json_out.is_file()
    assert markdown_out.is_file()

    payload = json.loads(json_out.read_text(encoding="utf-8"))
    assert payload["decision"]["recommendation"] == "defer"
    assert payload["decision"]["decision"] == "insufficient_data"
    assert "Decision: defer (insufficient_data)" in result.stderr


def test_runner_rejects_less_than_three_real_datasets(tmp_path: Path) -> None:
    script = _script_path()
    real_paths = _create_real_dataset_files(tmp_path, count=2)

    cmd = [
        sys.executable,
        str(script),
        "--iterations",
        "1",
        "--synthetic-work-dir",
        str(tmp_path / "synthetic"),
        "--size",
        "small=2",
        "--skip-polars",
        "--json-out",
        str(tmp_path / "result.json"),
        "--markdown-out",
        str(tmp_path / "result.md"),
    ]
    for idx, path in enumerate(real_paths, start=1):
        cmd.extend(["--real-dataset", f"real{idx}={path}"])

    result = subprocess.run(
        cmd,
        cwd=_repo_root(),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert "At least 3 real datasets are required" in result.stderr
