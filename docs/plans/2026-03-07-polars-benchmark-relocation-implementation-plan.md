# Polars Benchmark Relocation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Move the polars GFF benchmark helpers out of the runtime package into a repo-local `benchmarks/` package while keeping runner/test imports working and keeping the wheel clean.

**Architecture:** Perform a structural module move with no behavior changes. Rewrite/import tests first where needed, move the module into top-level `benchmarks/`, delete the old runtime copy, then verify both runner behavior and packaging boundaries.

**Tech Stack:** Python 3.13, pytest, Ruff, mypy, Hatchling, pathlib

---

### Task 1: Add Packaging Boundary Tests

**Files:**
- Modify: `tests/test_polars_gff_benchmark.py`
- Modify: `tests/test_polars_gff_benchmark_runner.py`
- Test: `tests/test_polars_gff_benchmark.py`
- Test: `tests/test_polars_gff_benchmark_runner.py`

**Step 1: Add a failing import-path test in `tests/test_polars_gff_benchmark.py`**

Add a test that imports `benchmarks.polars_gff_benchmark` instead of
`SpliceGrapher.formats.polars_gff_benchmark`.

**Step 2: Add a packaging-boundary assertion**

Add a small test that reads `pyproject.toml` and asserts the wheel target still
packages only `SpliceGrapher`.

**Step 3: Update the runner test import**

Change `tests/test_polars_gff_benchmark_runner.py` to import
`write_synthetic_gff` from `benchmarks.polars_gff_benchmark`.

**Step 4: Run the benchmark-focused test slice to confirm red**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
```

Expected: fail because `benchmarks.polars_gff_benchmark` does not exist yet.

**Step 5: Commit the failing test rewrite**

```bash
git add tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
git commit -m "test: point benchmark coverage at top-level package"
```

### Task 2: Move The Benchmark Module

**Files:**
- Create: `benchmarks/__init__.py`
- Create: `benchmarks/polars_gff_benchmark.py`
- Delete: `SpliceGrapher/formats/polars_gff_benchmark.py`

**Step 1: Create `benchmarks/__init__.py`**

Create an empty package marker or a minimal module docstring.

**Step 2: Move the module**

Copy the contents of `SpliceGrapher/formats/polars_gff_benchmark.py` into
`benchmarks/polars_gff_benchmark.py`.

**Step 3: Delete the old runtime module**

Remove `SpliceGrapher/formats/polars_gff_benchmark.py` entirely. Do not leave a
shim or re-export behind.

**Step 4: Run Ruff format/check on the moved module**

Run:

```bash
uv run ruff format benchmarks/polars_gff_benchmark.py
uv run ruff check benchmarks/polars_gff_benchmark.py
```

Expected: pass.

**Step 5: Commit the move**

```bash
git add benchmarks/__init__.py benchmarks/polars_gff_benchmark.py
git rm SpliceGrapher/formats/polars_gff_benchmark.py
git commit -m "refactor: move polars benchmark module out of runtime package"
```

### Task 3: Rewire Runner Imports

**Files:**
- Modify: `scripts/benchmarks/run_gff_loader_benchmarks.py`
- Test: `tests/test_polars_gff_benchmark_runner.py`

**Step 1: Update the runner import**

Replace:

```python
from SpliceGrapher.formats.polars_gff_benchmark import ...
```

With:

```python
from benchmarks.polars_gff_benchmark import ...
```

**Step 2: Grep for leftover runtime imports**

Run:

```bash
rg -n "SpliceGrapher\\.formats\\.polars_gff_benchmark|from benchmarks\\.polars_gff_benchmark" SpliceGrapher tests scripts
```

Expected: no runtime imports remain; only the new top-level package import is left.

**Step 3: Run the benchmark-focused test slice**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
```

Expected: pass.

**Step 4: Commit the import migration**

```bash
git add scripts/benchmarks/run_gff_loader_benchmarks.py tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
git commit -m "refactor: rewire benchmark imports to top-level package"
```

### Task 4: Final Verification

**Files:**
- Verify: `benchmarks/polars_gff_benchmark.py`
- Verify: `scripts/benchmarks/run_gff_loader_benchmarks.py`
- Verify: `tests/test_polars_gff_benchmark.py`
- Verify: `tests/test_polars_gff_benchmark_runner.py`
- Verify: `pyproject.toml`

**Step 1: Run final touched checks**

Run:

```bash
uv run ruff check benchmarks/polars_gff_benchmark.py scripts/benchmarks/run_gff_loader_benchmarks.py tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
uv run ruff format --check benchmarks/polars_gff_benchmark.py scripts/benchmarks/run_gff_loader_benchmarks.py tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_polars_gff_benchmark.py tests/test_polars_gff_benchmark_runner.py
uv build
```

Expected: all pass.

**Step 2: Run a broader regression slice if the touched checks are clean**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider
```

Expected: pass.

**Step 3: Commit any final cleanups**

```bash
git add benchmarks scripts tests
git commit -m "test: verify benchmark package relocation boundary"
```
