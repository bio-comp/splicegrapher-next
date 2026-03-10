# Shared Process Hard Cut Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the legacy shared process helper module with a modern snake_case module and remove the remaining Python 2 residue from the shared progress helper.

**Architecture:** Delete `SpliceGrapher/shared/process_utils.py` and replace it with `SpliceGrapher/shared/process.py`, updating all SGN imports and tests in the same PR. Keep `SpliceGrapher/shared/progress.py` in place, but remove Python 2 iterator residue and modernize any touched interfaces without adding compatibility aliases.

**Tech Stack:** Python 3.11+, `structlog`, `subprocess`, `pytest`, `ruff`, `mypy`

---

### Task 1: Pin the new shared process boundary with failing tests

**Files:**
- Modify: `tests/test_process_utils.py`
- Modify: `tests/test_core_shared_import_smoke.py`
- Modify: `tests/test_shared_utils_shim.py`

**Step 1: Write the failing tests**

Rename the process helper expectations to the new boundary:

```python
from SpliceGrapher.shared import process as process_module


def test_run_logged_command_signature_is_explicit() -> None:
    signature = inspect.signature(process_module.run_logged_command)
    assert all(param.kind != inspect.Parameter.VAR_KEYWORD for param in signature.parameters.values())
```

Update smoke tests to import `SpliceGrapher.shared.process` and stop asserting `process_utils` exports.

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_process_utils.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py`

Expected: FAIL because the new `shared.process` module and renamed APIs do not exist yet.

**Step 3: Commit the failing tests**

```bash
git add tests/test_process_utils.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py
git commit -m "test: pin shared process boundary"
```

### Task 2: Create `shared/process.py` and remove the legacy process module

**Files:**
- Delete: `SpliceGrapher/shared/process_utils.py`
- Create: `SpliceGrapher/shared/process.py`
- Modify: `tests/test_process_utils.py`

**Step 1: Implement the new process module**

Create `SpliceGrapher/shared/process.py` with:
- `get_attribute`
- `id_factory`
- `log_message`
- `run_command`
- `run_logged_command`
- `write_startup_message`

Requirements:
- no camelCase names
- no `**args: object`
- no `sys.stderr.write(...)`
- `run_logged_command` should preserve old `runCommand` behavior where tests pin it

**Step 2: Run the focused process test slice**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_process_utils.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py`

Expected: PASS or near-pass with only progress-related fallout left.

**Step 3: Commit**

```bash
git add SpliceGrapher/shared/process.py SpliceGrapher/shared/process_utils.py tests/test_process_utils.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py
git commit -m "refactor(shared): replace legacy process helpers"
```

### Task 3: Update SGN runtime imports to the new shared process boundary

**Files:**
- Modify any SGN runtime callers found by `rg` for `process_utils`
- Modify any tests still importing `process_utils`

**Step 1: Rewrite imports and names**

Examples:

```python
from SpliceGrapher.shared.process import get_attribute, run_logged_command
```

Remove all remaining `process_utils` imports and all uses of:
- `getAttribute`
- `logMessage`
- `runCommand`
- `writeStartupMessage`

**Step 2: Run the focused related tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_process_utils.py tests/test_alignment_io_process_utils_boundary.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py`

Expected: PASS

**Step 3: Commit**

```bash
git add SpliceGrapher tests
git commit -m "refactor(shared): update callers to shared process module"
```

### Task 4: Remove Python 2 residue from `shared/progress.py`

**Files:**
- Modify: `SpliceGrapher/shared/progress.py`
- Modify: `tests/test_progress.py`

**Step 1: Write or adjust the failing progress test if needed**

Pin the absence of the Python 2 alias:

```python
def test_random_list_iterator_has_no_python2_next_alias() -> None:
    iterator = progress_module.RandomListIterator([1, 2, 3], seed=1)
    assert not hasattr(iterator, "next")
```

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_progress.py`

Expected: FAIL because `next` still exists.

**Step 3: Implement the minimal cleanup**

Remove `next = __next__` and modernize any touched naming only if required by the tests.

**Step 4: Run test to verify it passes**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_progress.py`

Expected: PASS

**Step 5: Commit**

```bash
git add SpliceGrapher/shared/progress.py tests/test_progress.py
git commit -m "refactor(shared): remove progress python2 residue"
```

### Task 5: Run focused shared-helper quality gates

**Files:**
- Verify only; no expected new files

**Step 1: Run Ruff**

Run: `uv run ruff check SpliceGrapher/shared/process.py SpliceGrapher/shared/progress.py tests/test_process_utils.py tests/test_progress.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py --fix`

Expected: PASS

**Step 2: Run formatting**

Run: `uv run ruff format SpliceGrapher/shared/process.py SpliceGrapher/shared/progress.py tests/test_process_utils.py tests/test_progress.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py`

Expected: PASS

**Step 3: Run MyPy**

Run: `uv run mypy SpliceGrapher/shared/process.py SpliceGrapher/shared/progress.py tests/test_process_utils.py tests/test_progress.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py`

Expected: PASS

**Step 4: Run focused tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_process_utils.py tests/test_progress.py tests/test_core_shared_import_smoke.py tests/test_shared_utils_shim.py tests/test_alignment_io_process_utils_boundary.py`

Expected: PASS

### Task 6: Run full repository verification

**Files:**
- Verify only; no expected new files

**Step 1: Run full Ruff**

Run: `uv run ruff check . --fix && uv run ruff format .`

Expected: PASS

**Step 2: Run full pytest**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`

Expected: PASS

**Step 3: Run full MyPy**

Run: `uv run mypy SpliceGrapher tests`

Expected: PASS

**Step 4: Run build**

Run: `uv build`

Expected: PASS

**Step 5: Final commit if verification changed formatting**

```bash
git add SpliceGrapher/shared/process.py SpliceGrapher/shared/progress.py tests
git commit -m "test: finalize shared helper hard cut"
```
