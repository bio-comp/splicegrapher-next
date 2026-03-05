"""Clean-invariant ratchet guard for SGN."""

from __future__ import annotations

import subprocess
import sys
import tomllib
from pathlib import Path, PurePosixPath

BASELINE_PER_FILE_IGNORES: dict[str, set[str]] = {
    "SpliceGrapher/formats/fasta.py": {"E", "F", "W"},
}

# Real script entry points can be added here explicitly as needed.
EXECUTABLE_PYTHON_ALLOWLIST: set[str] = set()

# New per-file ignores are only allowed for explicit non-production carve-outs.
ALLOWED_NEW_IGNORE_TARGETS: tuple[str, ...] = (
    "SpliceGrapher/shims/**/*.py",
    "tests/**/*.py",
    "docs/notebooks/*.ipynb",
)
ENUM_CONTROL_FLOW_PATHS: tuple[str, ...] = (
    "SpliceGrapher/SpliceGraph.py",
    "SpliceGrapher/formats/GeneModel.py",
)
MAGIC_STRING_CONTROL_FLOW_LITERALS: tuple[str, ...] = (
    "gene",
    "mrna",
    "exon",
    "cds",
    "predicted_gene",
    "pseudogene",
    "pseudogenic_transcript",
    "pseudogenic_exon",
    "parent",
    "child",
    "known",
    "predicted",
    "unresolved",
)
OVERLAP_CONTROL_FLOW_PATHS: tuple[str, ...] = (
    "SpliceGrapher/SpliceGraph.py",
    "SpliceGrapher/formats/GeneModel.py",
)
# Ratchet baseline: existing manual overlap logic is tolerated, but new instances
# in protected modules must route through interval helper abstractions.
BASELINE_MANUAL_OVERLAP_LINES: dict[str, set[str]] = {
    "SpliceGrapher/SpliceGraph.py": {
        "if e.minpos < n.minpos and n.maxpos < e.maxpos:",
        "return self.maxpos > o.minpos and self.minpos < o.maxpos",
    },
    "SpliceGrapher/formats/GeneModel.py": {
        "return strand == self.strand and self.minpos <= pos <= self.maxpos",
        "if g.maxpos < minpos or g.minpos > maxpos:",
    },
}


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


def _load_sources(repo_root: Path, relative_paths: tuple[str, ...]) -> dict[str, str]:
    source_by_path: dict[str, str] = {}
    for rel_path in relative_paths:
        file_path = repo_root / rel_path
        if file_path.exists():
            source_by_path[rel_path] = file_path.read_text(encoding="utf-8")
    return source_by_path


def find_magic_string_control_flow(
    *,
    source_by_path: dict[str, str],
    protected_paths: tuple[str, ...],
) -> list[str]:
    """Flag direct literal-string branching in enum-protected modules."""
    violations: list[str] = []
    for rel_path in protected_paths:
        source = source_by_path.get(rel_path, "")
        for line_no, line in enumerate(source.splitlines(), start=1):
            stripped = line.strip()
            if stripped.startswith("#"):
                continue
            if not (stripped.startswith("if ") or stripped.startswith("elif ")):
                continue
            for literal in MAGIC_STRING_CONTROL_FLOW_LITERALS:
                if f'"{literal}"' in stripped or f"'{literal}'" in stripped:
                    violations.append(
                        f"{rel_path}:{line_no}: magic-string control-flow literal {literal!r}"
                    )
                    break
    return violations


def _looks_like_manual_overlap_line(stripped: str) -> bool:
    if stripped.startswith("#"):
        return False
    if not ("if " in stripped or "elif " in stripped or stripped.startswith("return ")):
        return False
    if "minpos" not in stripped or "maxpos" not in stripped:
        return False
    if " and " not in stripped and " or " not in stripped:
        return False
    if not any(token in stripped for token in ("<", ">", "<=", ">=")):
        return False
    if "intervals_overlap(" in stripped or "interval_contains(" in stripped:
        return False
    return True


def find_manual_overlap_lines(
    *,
    source_by_path: dict[str, str],
    protected_paths: tuple[str, ...],
) -> dict[str, set[str]]:
    """Collect normalized manual overlap control-flow lines per protected path."""
    result: dict[str, set[str]] = {}
    for rel_path in protected_paths:
        source = source_by_path.get(rel_path, "")
        for line in source.splitlines():
            stripped = line.strip()
            if not _looks_like_manual_overlap_line(stripped):
                continue
            result.setdefault(rel_path, set()).add(stripped)
    return result


def detect_manual_overlap_growth(
    *,
    current: dict[str, set[str]],
    baseline: dict[str, set[str]],
) -> list[str]:
    """Return violations where manual overlap logic expands in protected modules."""
    violations: list[str] = []
    for path in sorted(current):
        known = baseline.get(path, set())
        unexpected = sorted(current[path] - known)
        for line in unexpected:
            violations.append(f"{path}: unexpected manual-overlap control flow -> {line}")
    return violations


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
    source_by_path = _load_sources(repo_root, ENUM_CONTROL_FLOW_PATHS)
    magic_string_violations = find_magic_string_control_flow(
        source_by_path=source_by_path,
        protected_paths=ENUM_CONTROL_FLOW_PATHS,
    )
    overlap_sources = _load_sources(repo_root, OVERLAP_CONTROL_FLOW_PATHS)
    current_manual_overlap = find_manual_overlap_lines(
        source_by_path=overlap_sources,
        protected_paths=OVERLAP_CONTROL_FLOW_PATHS,
    )
    manual_overlap_violations = detect_manual_overlap_growth(
        current=current_manual_overlap,
        baseline=BASELINE_MANUAL_OVERLAP_LINES,
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

    if magic_string_violations:
        sys.stdout.write(
            "Clean-invariant failure: magic-string control flow found in enum-protected modules.\n"
        )
        for violation in magic_string_violations:
            sys.stdout.write(f"- {violation}\n")

    if manual_overlap_violations:
        sys.stdout.write(
            "Clean-invariant failure: new manual overlap control flow found in protected modules.\n"
        )
        for violation in manual_overlap_violations:
            sys.stdout.write(f"- {violation}\n")

    if (
        ignore_violations
        or executable_violations
        or magic_string_violations
        or manual_overlap_violations
    ):
        return 1

    sys.stdout.write("Clean-invariant checks passed.\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
