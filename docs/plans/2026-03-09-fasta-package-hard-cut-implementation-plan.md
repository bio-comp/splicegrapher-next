# FASTA Package Hard Cut Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the legacy flat FASTA module with a modern `SpliceGrapher.formats.fasta` package, delete Python 2 residue, and expose a snake_case-only public API.

**Architecture:** Convert `SpliceGrapher/formats/fasta.py` into a package with explicit module boundaries for records, readers, and operations. Keep the top-level `SpliceGrapher.formats.fasta` namespace stable while deleting legacy names and updating SGN tests to the new modern API.

**Tech Stack:** Python 3.11+, `structlog`, `pyfaidx`, `pytest`, `ruff`, `mypy`

---

### Task 1: Add the package layout guard and new public API expectations

**Files:**
- Modify: `tests/test_fasta_io.py`
- Create: `tests/test_fasta_module_layout.py`

**Step 1: Write the failing tests**

Add a new layout test that imports `SpliceGrapher.formats.fasta` and asserts:

```python
from SpliceGrapher.formats import fasta


def test_fasta_is_package_with_modern_exports() -> None:
    assert hasattr(fasta, "FastaIterator")
    assert hasattr(fasta, "FastaSlice")
    assert hasattr(fasta, "truncate_sequences")
    assert not hasattr(fasta, "fasta_itr")
    assert not hasattr(fasta, "fasta_slice")
    assert not hasattr(fasta, "truncateSequences")
```

Update `tests/test_fasta_io.py` to use the new names:

```python
records = list(fasta.FastaIterator(fasta_path))
```

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: FAIL because the new package/classes do not exist yet.

**Step 3: Commit the failing test**

```bash
git add tests/test_fasta_io.py tests/test_fasta_module_layout.py
git commit -m "test: pin FASTA package boundary"
```

### Task 2: Create the package boundary and move record types

**Files:**
- Delete: `SpliceGrapher/formats/fasta.py`
- Create: `SpliceGrapher/formats/fasta/__init__.py`
- Create: `SpliceGrapher/formats/fasta/records.py`

**Step 1: Write minimal record modules**

Create `records.py` with:

```python
class MalformedInput(Exception):
    ...


@dataclass(frozen=True, slots=True)
class FastaRecord:
    header: str
    sequence: str
```

Create `__init__.py` that re-exports the modern public boundary.

**Step 2: Run the focused tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_module_layout.py`

Expected: still FAIL because reader/operation classes are not implemented yet, but import/package structure should resolve.

**Step 3: Commit**

```bash
git add SpliceGrapher/formats/fasta/__init__.py SpliceGrapher/formats/fasta/records.py
git commit -m "refactor(formats): add FASTA package boundary"
```

### Task 3: Move reader logic into a typed Python 3 reader module

**Files:**
- Create: `SpliceGrapher/formats/fasta/readers.py`
- Modify: `SpliceGrapher/formats/fasta/__init__.py`

**Step 1: Implement the low-level readers**

Move and modernize:

- `_to_text`
- `_fasta_itr_from_file`
- `_fasta_itr_from_pyfaidx`
- `_fasta_itr_from_name`
- `_fasta_itr`
- `fasta_get_by_name`
- `FastaIterator`
- `FastaSlice`

Modernization requirements:

- no `(object)` bases
- no `next = __next__`
- typed constructor signatures
- snake_case parameter names
- no legacy alias exports

**Step 2: Run the focused FASTA tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: partial PASS; operation-related failures may remain.

**Step 3: Commit**

```bash
git add SpliceGrapher/formats/fasta/readers.py SpliceGrapher/formats/fasta/__init__.py
git commit -m "refactor(formats): modernize FASTA readers"
```

### Task 4: Move operation helpers into a modern operations module

**Files:**
- Create: `SpliceGrapher/formats/fasta/operations.py`
- Modify: `SpliceGrapher/formats/fasta/__init__.py`
- Modify: `tests/test_fasta_io.py`

**Step 1: Implement operations**

Move and modernize:

- `get_sequence`
- `fasta_count`
- `fasta_split`
- `FastaRandomizer`
- `truncate_sequences`

Requirements:

- rename `truncateSequences` to `truncate_sequences`
- rename all camelCase parameters/local variables to snake_case
- replace `%` interpolation with f-strings
- replace direct `sys.stderr.write(...)` with `structlog`
- avoid `typing.Any`

**Step 2: Run the FASTA test slice**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: PASS

**Step 3: Commit**

```bash
git add SpliceGrapher/formats/fasta/operations.py SpliceGrapher/formats/fasta/__init__.py tests/test_fasta_io.py tests/test_fasta_module_layout.py
git commit -m "refactor(formats): modernize FASTA operations"
```

### Task 5: Run quality gates on the touched FASTA slice

**Files:**
- Verify only; no expected new files

**Step 1: Run Ruff**

Run: `uv run ruff check SpliceGrapher/formats/fasta tests/test_fasta_io.py tests/test_fasta_module_layout.py --fix`

Expected: PASS

**Step 2: Run formatting**

Run: `uv run ruff format SpliceGrapher/formats/fasta tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: PASS

**Step 3: Run MyPy**

Run: `uv run mypy SpliceGrapher/formats/fasta tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: PASS

**Step 4: Run focused tests**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_io.py tests/test_fasta_module_layout.py`

Expected: PASS

### Task 6: Run the repository-level verification gates

**Files:**
- Verify only; no expected new files

**Step 1: Run full pytest**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`

Expected: PASS

**Step 2: Run build**

Run: `uv build`

Expected: PASS

**Step 3: Final commit if verification changed formatting**

```bash
git add SpliceGrapher/formats/fasta tests/test_fasta_io.py tests/test_fasta_module_layout.py
git commit -m "test: finalize FASTA package hard cut"
```
