# Annotation IO Hard-Cut Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the loose `annotation_io.load_gene_models(..., **args)` boundary with an explicit typed API while keeping the module in place and preserving behavior.

**Architecture:** Keep `SpliceGrapher/formats/annotation_io.py` as a single module, but reorganize it into clearer internal sections and make `load_gene_models(...)` the only explicit public boundary. The hard cut is at the function signature and typed option flow, not at the file layout.

**Tech Stack:** Python 3.12, gffutils, pathlib, pytest, Ruff, MyPy

---

### Task 1: Pin the new `load_gene_models` boundary in tests

**Files:**
- Modify: `tests/test_annotation_io.py`
- Test: `tests/test_annotation_io.py`

**Step 1: Write the failing test**

Add a test that calls `load_gene_models(..., outdir=..., cache_dir=...)` and a test that asserts an unknown keyword raises `TypeError`.

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_annotation_io.py`
Expected: FAIL because the old permissive `**args` boundary still exists and the new explicit contract is not pinned.

**Step 3: Write minimal test updates**

Update `tests/test_annotation_io.py` to assert:
- explicit named arguments are accepted
- unknown keyword arguments are rejected

**Step 4: Run test to verify it passes or fails for the right reason**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_annotation_io.py`
Expected: FAIL only because implementation still needs the explicit signature.

**Step 5: Commit**

```bash
git add tests/test_annotation_io.py
git commit -m "test: pin annotation_io boundary"
```

### Task 2: Hard-cut `load_gene_models(...)` to an explicit typed signature

**Files:**
- Modify: `SpliceGrapher/formats/annotation_io.py`
- Test: `tests/test_annotation_io.py`

**Step 1: Replace the loose function signature**

Change:
- `load_gene_models(path: str, **args: object) -> GeneModel`

To:
- `load_gene_models(path: str | Path, *, outdir: str | Path | None = None, cache_dir: str | Path | None = None) -> GeneModel`

**Step 2: Remove `cast(...)`-driven option extraction**

Replace the internal `args.get(...)` flow with direct named parameters and `Path(...)` normalization where needed.

**Step 3: Keep behavior stable**

Preserve:
- path existence checks
- optional DB cache creation under `cache_dir`
- optional intron BED cache writing under `outdir`

**Step 4: Run focused tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_annotation_io.py tests/test_integration_simple.py`
Expected: PASS

**Step 5: Commit**

```bash
git add SpliceGrapher/formats/annotation_io.py tests/test_annotation_io.py
git commit -m "refactor(formats): hard-cut annotation_io boundary"
```

### Task 3: Tighten internal helper typing and layout in-place

**Files:**
- Modify: `SpliceGrapher/formats/annotation_io.py`
- Test: `tests/test_annotation_io.py`

**Step 1: Reorganize helper sections**

Keep the file in place, but group it into clear blocks:
- feature normalization / attribute helpers
- DB/cache helpers
- gene-model assembly helpers
- intron cache helpers
- public loader

**Step 2: Reduce loose typing**

Prefer concrete types over cast-heavy flows where possible, especially around path handling and normalized gene records.

**Step 3: Run static checks on touched files**

Run:
- `uv run ruff check SpliceGrapher/formats/annotation_io.py tests/test_annotation_io.py tests/test_integration_simple.py`
- `uv run mypy SpliceGrapher/formats/annotation_io.py tests/test_annotation_io.py tests/test_integration_simple.py`

Expected: PASS

**Step 4: Commit**

```bash
git add SpliceGrapher/formats/annotation_io.py
git commit -m "refactor(formats): tighten annotation_io internals"
```

### Task 4: Run full verification and prepare the PR

**Files:**
- Verify repo state only

**Step 1: Run repo-wide gates**

Run:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

Expected: all pass

**Step 2: Inspect git status**

Confirm only intended code/test/docs files are staged. Keep `AGENTS.md` and `LOCAL_TODO.md` unstaged.

**Step 3: Open the PR**

Describe:
- the explicit `load_gene_models(...)` contract
- removal of permissive `**args`
- migration impact for downstream callers

**Step 4: Label the PR**

Apply:
- `area:formats`
- `type:refactor`
- `type:test` if tests changed
