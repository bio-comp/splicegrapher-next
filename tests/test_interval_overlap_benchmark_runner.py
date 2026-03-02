from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _script_path() -> Path:
    return _repo_root() / "scripts" / "benchmarks" / "run_interval_overlap_benchmarks.py"


def test_interval_overlap_benchmark_runner_writes_outputs(tmp_path: Path) -> None:
    script = _script_path()
    json_out = tmp_path / "result.json"
    markdown_out = tmp_path / "result.md"

    cmd = [
        sys.executable,
        str(script),
        "--interval-count",
        "200",
        "--query-count",
        "200",
        "--iterations",
        "1",
        "--assert-threshold",
        "--max-slowdown-ratio",
        "5.0",
        "--json-out",
        str(json_out),
        "--markdown-out",
        str(markdown_out),
    ]

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
    assert payload["parity_ok"] is True
    assert payload["indexed_slowdown_ratio"] > 0.0
    assert "Threshold check: pass" in result.stderr


def test_interval_overlap_benchmark_runner_fails_on_strict_threshold(tmp_path: Path) -> None:
    script = _script_path()

    cmd = [
        sys.executable,
        str(script),
        "--interval-count",
        "200",
        "--query-count",
        "200",
        "--iterations",
        "1",
        "--assert-threshold",
        "--max-slowdown-ratio",
        "0.01",
        "--json-out",
        str(tmp_path / "result.json"),
        "--markdown-out",
        str(tmp_path / "result.md"),
    ]

    result = subprocess.run(
        cmd,
        cwd=_repo_root(),
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode != 0
    assert "Threshold check: fail" in result.stderr
