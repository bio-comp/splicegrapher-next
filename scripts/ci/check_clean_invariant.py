"""Clean-invariant ratchet guard for SGN."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path, PurePosixPath

try:
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python < 3.11 fallback
    import tomli as tomllib


BASELINE_PER_FILE_IGNORES: dict[str, set[str]] = {
    "SpliceGrapher/SpliceGraph.py": {"E", "F"},
    "SpliceGrapher/formats/GeneModel.py": {"E", "F", "W"},
    "SpliceGrapher/formats/alignment_io.py": {"E", "F", "W"},
    "SpliceGrapher/formats/fasta.py": {"E", "F", "W"},
    "SpliceGrapher/shared/ShortRead.py": {"E", "F", "W"},
}

# Real script entry points can be added here explicitly as needed.
EXECUTABLE_PYTHON_ALLOWLIST: set[str] = set()

# New per-file ignores are only allowed for explicit non-production carve-outs.
ALLOWED_NEW_IGNORE_TARGETS: tuple[str, ...] = (
    "SpliceGrapher/shims/**/*.py",
    "tests/**/*.py",
    "docs/notebooks/*.ipynb",
)


def _is_allowed_new_ignore_target(path: str) -> bool:
    pure_path = PurePosixPath(path)
    return any(pure_path.match(pattern) for pattern in ALLOWED_NEW_IGNORE_TARGETS)


def load_per_file_ignores(pyproject_path: Path) -> dict[str, set[str]]:
    """Load Ruff per-file-ignore configuration from pyproject."""
    with pyproject_path.open("rb") as handle:
        pyproject_data = tomllib.load(handle)
    per_file = (
        pyproject_data.get("tool", {}).get("ruff", {}).get("lint", {}).get("per-file-ignores", {})
    )
    return {str(path): set(values) for path, values in per_file.items()}


def detect_ignore_growth(
    *,
    current: dict[str, set[str]],
    baseline: dict[str, set[str]],
) -> list[str]:
    """Return ratchet violations where ignore debt grows."""
    violations: list[str] = []

    for path in sorted(current):
        if path not in baseline:
            if _is_allowed_new_ignore_target(path):
                continue
            violations.append(
                "New per-file ignore target introduced outside allowed carve-outs: "
                f"{path} -> {sorted(current[path])}"
            )
            continue

        extra_codes = sorted(current[path] - baseline[path])
        if extra_codes:
            violations.append(
                f"Per-file ignores widened for {path}: unexpected codes {extra_codes}"
            )

    return violations


def find_unexpected_executable_python_paths(
    *,
    ls_files_lines: list[str],
    executable_allowlist: set[str],
) -> list[str]:
    """Find executable .py paths not explicitly allowlisted."""
    unexpected: list[str] = []

    for line in ls_files_lines:
        if "\t" not in line:
            continue
        meta, path = line.split("\t", 1)
        mode = meta.split(" ", 1)[0]

        if mode != "100755":
            continue
        if not path.endswith(".py"):
            continue
        if path in executable_allowlist:
            continue

        unexpected.append(path)

    return sorted(unexpected)


def _tracked_python_mode_lines(repo_root: Path) -> list[str]:
    completed = subprocess.run(
        ["git", "-C", str(repo_root), "ls-files", "-s", "*.py"],
        check=True,
        capture_output=True,
        text=True,
    )
    return [line for line in completed.stdout.splitlines() if line.strip()]


def main() -> int:
    repo_root = Path(__file__).resolve().parents[2]
    current_ignores = load_per_file_ignores(repo_root / "pyproject.toml")
    ignore_violations = detect_ignore_growth(
        current=current_ignores,
        baseline=BASELINE_PER_FILE_IGNORES,
    )

    ls_mode_lines = _tracked_python_mode_lines(repo_root)
    executable_violations = find_unexpected_executable_python_paths(
        ls_files_lines=ls_mode_lines,
        executable_allowlist=EXECUTABLE_PYTHON_ALLOWLIST,
    )

    if ignore_violations:
        sys.stdout.write("Clean-invariant ratchet failure: per-file-ignore debt grew.\n")
        for violation in ignore_violations:
            sys.stdout.write(f"- {violation}\n")

    if executable_violations:
        sys.stdout.write(
            "Clean-invariant failure: unexpected executable Python modules detected.\n"
        )
        for path in executable_violations:
            sys.stdout.write(f"- {path}\n")

    if ignore_violations or executable_violations:
        return 1

    sys.stdout.write("Clean-invariant checks passed.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
