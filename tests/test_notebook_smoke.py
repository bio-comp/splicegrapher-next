"""Notebook execution smoke tests for SGN tutorial docs."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import nbformat
import pytest
from nbclient import NotebookClient

REPO_ROOT = Path(__file__).resolve().parents[1]
NOTEBOOK_PATHS = [
    REPO_ROOT / "docs" / "notebooks" / "01_legacy_tutorial_foundation.ipynb",
    REPO_ROOT / "docs" / "notebooks" / "02_modern_dataset_track.ipynb",
]


def _ensure_python3_kernel() -> None:
    subprocess.run(
        [
            sys.executable,
            "-m",
            "ipykernel",
            "install",
            "--user",
            "--name",
            "python3",
            "--display-name",
            "Python 3",
        ],
        check=True,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
    )


@pytest.fixture(scope="session", autouse=True)
def _register_kernel_once() -> None:
    _ensure_python3_kernel()


def _execute_notebook(
    path: Path,
    execution_cwd: Path = REPO_ROOT,
) -> nbformat.NotebookNode:
    notebook = nbformat.read(path, as_version=4)
    client = NotebookClient(
        notebook,
        timeout=120,
        kernel_name="python3",
        resources={"metadata": {"path": str(execution_cwd)}},
    )
    return client.execute()


def _text_outputs(notebook: nbformat.NotebookNode) -> str:
    chunks: list[str] = []
    for cell in notebook.cells:
        if cell.get("cell_type") != "code":
            continue
        for output in cell.get("outputs", []):
            if output.get("output_type") == "stream":
                chunks.append(str(output.get("text", "")))
                continue
            data = output.get("data", {})
            if "text/plain" in data:
                chunks.append(str(data["text/plain"]))
    return "\n".join(chunks)


@pytest.mark.integration
@pytest.mark.parametrize("notebook_path", NOTEBOOK_PATHS, ids=lambda p: p.name)
def test_tutorial_notebook_executes(notebook_path: Path) -> None:
    executed = _execute_notebook(notebook_path)
    output_blob = _text_outputs(executed)
    assert "ath" in output_blob
    assert "hsa" in output_blob


@pytest.mark.integration
def test_notebook_executes_from_notebook_directory_cwd() -> None:
    executed = _execute_notebook(
        NOTEBOOK_PATHS[0],
        execution_cwd=REPO_ROOT / "docs" / "notebooks",
    )
    output_blob = _text_outputs(executed)
    assert "PRJNA755474" in output_blob
