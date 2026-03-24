# Alignment Collect Lean Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** make `SpliceGrapher/formats/alignment_io/collect.py` a pure pysam-backed collection engine without changing the public `alignment_io` API.

**Architecture:** keep the existing `alignment_io` package shape. Move source classification and depths-file fallback bridging into `api.py`, and leave `collect.py` with only pysam traversal, junction construction, and depth-array helpers. Do not add new internal files.

**Tech Stack:** Python 3.13, `pysam`, `numpy`, `pytest`, `ruff`, `mypy`, `uv`

---

### Task 1: Pin the lean boundary in tests

**Files:**
- Modify: `tests/test_alignment_io_process_utils_boundary.py`
- Modify: `tests/test_alignment_io_constants.py`
- Test: `tests/test_alignment_io_process_utils_boundary.py`
- Test: `tests/test_alignment_io_constants.py`

**Step 1: Write the failing test**

Update the boundary tests so they expect:
- `alignment_collect` still exposes `_collect_pysam_data`
- `alignment_collect` no longer exposes `_is_depths_source`
- `alignment_collect` no longer exposes `_collect_depths_source_data`
- `alignment_api` now owns `_is_depths_source`
- focused fallback tests patch `alignment_api._is_depths_source` and `alignment_api._collect_depths_source_data`

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py`

Expected: FAIL because the helpers still live in `collect.py`.

**Step 3: Write minimal implementation**

Only land the red expectations in the tests.

**Step 4: Run test to verify it fails cleanly**

Run the same command again.

Expected: FAIL for missing or stale helper ownership.

**Step 5: Commit**

```bash
git add tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py
git commit -m "test: pin alignment collect lean boundary"
```

### Task 2: Move the non-pysam helpers into `api.py`

**Files:**
- Modify: `SpliceGrapher/formats/alignment_io/api.py`
- Modify: `SpliceGrapher/formats/alignment_io/collect.py`
- Test: `tests/test_alignment_io_process_utils_boundary.py`
- Test: `tests/test_alignment_io_constants.py`

**Step 1: Write the failing test**

Use the red boundary tests from Task 1.

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py`

Expected: FAIL because `api.py` does not yet own the fallback helpers.

**Step 3: Write minimal implementation**

In `api.py`:
- add private `_is_depths_source(...)`
- add private `_collect_depths_source_data(...)`
- update `read_alignment_depths(...)`, `read_alignment_junctions(...)`, and `collect_alignment_data(...)` to call the new local helpers

In `collect.py`:
- remove `_is_depths_source(...)`
- remove `_collect_depths_source_data(...)`
- remove `depth_io` imports used only by those helpers
- keep only pysam-backed collection and its tight internal helpers
- shrink `__all__` accordingly

**Step 4: Run test to verify it passes**

Run the same command again.

Expected: PASS.

**Step 5: Commit**

```bash
git add SpliceGrapher/formats/alignment_io/api.py SpliceGrapher/formats/alignment_io/collect.py tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py
git commit -m "refactor(formats): hard-clean alignment collect boundary"
```

### Task 3: Prove parity and full repo stability

**Files:**
- Test: `tests/test_alignment_io_parity.py`
- Test: `tests/test_shortread_compat.py`
- Test: repo-wide verification commands

**Step 1: Run focused alignment tests**

Run:
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py tests/test_alignment_io_parity.py tests/test_shortread_compat.py`

Expected: PASS.

**Step 2: Run full verification**

Run:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

Expected: PASS.

**Step 3: Commit any final cleanup**

If verification required final import or formatting cleanup:

```bash
git add SpliceGrapher/formats/alignment_io/api.py SpliceGrapher/formats/alignment_io/collect.py tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_constants.py
git commit -m "refactor(formats): finalize alignment collect cleanup"
```
