# GeneModel Internal Split Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Split `SpliceGrapher/formats/gene_model/model.py` into smaller internal modules for query logic and loading/index maintenance while keeping `SpliceGrapher.formats.gene_model.GeneModel` behavior and imports stable.

**Architecture:** `GeneModel` remains the public class in `model.py`, but its larger internal responsibilities move into `queries.py` and `loading.py`. The class becomes thinner and delegates to focused helpers without changing parser, writer, or package semantics.

**Tech Stack:** Python 3.11+, dataclasses, `pytest`, `ruff`, `mypy`, Hatchling/`uv`

---

### Task 1: Pin the New Internal GeneModel Layout in Tests

**Files:**
- Modify: `tests/test_gene_model_package_layout.py`
- Modify: `tests/test_gene_model.py`
- Test: `tests/test_gene_model_package_layout.py`
- Test: `tests/test_gene_model.py`

**Step 1: Write the failing test**

Add assertions that expect:
- `SpliceGrapher.formats.gene_model.queries` to exist
- `SpliceGrapher.formats.gene_model.loading` to exist
- `SpliceGrapher.formats.gene_model.GeneModel` to remain the public owning class
- any signature/monkeypatch tests still resolve against the public class, not direct internals

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_package_layout.py tests/test_gene_model.py`

Expected: FAIL because the new internal modules do not exist yet.

**Step 3: Write minimal implementation**

Only land the red expectations in the tests.

**Step 4: Run test to verify it fails cleanly**

Run the same command again.

Expected: FAIL for missing-module or stale-layout reasons.

**Step 5: Commit**

```bash
git add tests/test_gene_model_package_layout.py tests/test_gene_model.py
git commit -m "test: pin gene model internal module boundary"
```

### Task 2: Extract Query Helpers into `queries.py`

**Files:**
- Create: `SpliceGrapher/formats/gene_model/queries.py`
- Modify: `SpliceGrapher/formats/gene_model/model.py`
- Test: `tests/test_gene_model.py`
- Test: `tests/test_gene_model_package_layout.py`

**Step 1: Write the failing test**

Use the red layout test from Task 1 and any targeted query-behavior tests in `tests/test_gene_model.py` that exercise:
- `get_gene`
- `get_gene_by_name`
- `get_genes_in_range`
- acceptor/donor aggregation
- `isoform_dict`

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model.py -k "get_gene or get_genes_in_range or isoform_dict or acceptor or donor" tests/test_gene_model_package_layout.py`

Expected: FAIL because `queries.py` does not exist and `model.py` has not delegated yet.

**Step 3: Write minimal implementation**

Create `queries.py` and move query-only helpers there. Keep `GeneModel` methods as thin delegators or wrappers so the public API does not change.

**Step 4: Run test to verify it passes**

Run the same command again.

Expected: PASS.

**Step 5: Commit**

```bash
git add SpliceGrapher/formats/gene_model/queries.py SpliceGrapher/formats/gene_model/model.py tests/test_gene_model.py tests/test_gene_model_package_layout.py
git commit -m "refactor(formats): extract gene model query helpers"
```

### Task 3: Extract Loading and Index Helpers into `loading.py`

**Files:**
- Create: `SpliceGrapher/formats/gene_model/loading.py`
- Modify: `SpliceGrapher/formats/gene_model/model.py`
- Test: `tests/test_gene_model.py`
- Test: `tests/test_gene_model_package_layout.py`

**Step 1: Write the failing test**

Keep the layout test red and target behavior around:
- `GeneModel.from_gff`
- `load_gene_model`
- `make_sorted_model`
- constructor/load failure handling for missing or empty inputs

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model.py -k "from_gff or load_gene_model or make_sorted_model or No gene models found" tests/test_gene_model_package_layout.py`

Expected: FAIL because `loading.py` does not exist and the delegation split is incomplete.

**Step 3: Write minimal implementation**

Create `loading.py` and move loading/index helpers there. Keep `GeneModel` as the owning state object and delegate loading/index operations through small methods.

**Step 4: Run test to verify it passes**

Run the same command again.

Expected: PASS.

**Step 5: Commit**

```bash
git add SpliceGrapher/formats/gene_model/loading.py SpliceGrapher/formats/gene_model/model.py tests/test_gene_model.py tests/test_gene_model_package_layout.py
git commit -m "refactor(formats): extract gene model loading helpers"
```

### Task 4: Tighten Package Layout, Verify Full Suite, and Update Present-State Docs

**Files:**
- Modify: `SpliceGrapher/formats/gene_model/__init__.py`
- Modify: `LOCAL_TODO.md` (local only, not committed)
- Test: repo-wide verification commands

**Step 1: Write the failing test**

No new behavior test is required beyond the layout and behavior checks already added, but use the package-layout test to pin any final internal export expectations.

**Step 2: Run verification before final cleanup**

Run:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

Expected: PASS or expose any missed import/delegation fallout.

**Step 3: Write minimal implementation**

Make any final import cleanups in `__init__.py` needed to keep the public package stable without turning it into another monolith.

**Step 4: Run verification to verify it still passes**

Run the full command set again if final cleanup changed tracked files.

**Step 5: Commit**

```bash
git add SpliceGrapher/formats/gene_model/__init__.py
git commit -m "refactor(formats): finalize gene model internal split"
```
