# Parser Facade Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Tighten `SpliceGrapher.formats.parsers` into an explicit package facade and modernize `SpliceGraphParser` by deleting `loadFromFile` and its manual iterator bookkeeping.

**Architecture:** Keep eager parsing in `SpliceGraphParser.__init__`, but rename the loader to `load_from_file`, remove `__next__`/index state, and expose `load_gene_model_records` plus `SpliceGraphParser` directly from `SpliceGrapher.formats.parsers`. Update tests and any fallout in the same branch.

**Tech Stack:** Python 3.12, `uv`, `pytest`, `ruff`, `mypy`

---

### Task 1: Write the failing parser-boundary tests

**Files:**
- Modify: `tests/test_gene_model_gff_module_layout.py`
- Modify: `tests/test_gene_model_gff_parser.py`
- Modify: `tests/test_splice_graph.py`

**Step 1: Tighten the parser boundary tests**

Update tests so they assert:
- `from SpliceGrapher.formats.parsers import load_gene_model_records, SpliceGraphParser` works
- `SpliceGraphParser` remains iterable via `iter(parser)` / `list(parser)` semantics
- if pinned explicitly, `load_from_file` exists and `loadFromFile` does not

**Step 2: Run the targeted tests to verify they fail**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_splice_graph.py
```

Expected: fail because `parsers.__init__` does not yet export the new boundary and `SpliceGraphParser` still exposes the legacy loader/iterator shape.

---

### Task 2: Implement the facade and iterator cleanup

**Files:**
- Modify: `SpliceGrapher/formats/parsers/__init__.py`
- Modify: `SpliceGrapher/formats/parsers/splice_graph.py`

**Step 1: Replace the package facade**

`SpliceGrapher/formats/parsers/__init__.py` should expose:

```python
from .gene_model_gff import load_gene_model_records
from .splice_graph import SpliceGraphParser
```

with matching `__all__`.

**Step 2: Modernize `SpliceGraphParser`**

In `SpliceGrapher/formats/parsers/splice_graph.py`:
- rename `loadFromFile` to `load_from_file`
- update `__init__` to call `self.load_from_file()`
- remove `self._graph_keys`
- remove `self.graphId`
- remove `__next__`
- change `__iter__` to return `Iterator[SpliceGraph]` over `self.graphDict.values()`

**Step 3: Fix any direct fallout from the renamed loader**

Search for `loadFromFile` in the repo and update or delete the callers/tests in this same task.

**Step 4: Run the targeted tests again**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_splice_graph.py
```

Expected: pass.

---

### Task 3: Run targeted lint/type/gene-model verification

**Files:**
- Modify only if fallout requires it

**Step 1: Run targeted lint and type checks**

Run:

```bash
uv run ruff check SpliceGrapher/formats/parsers/__init__.py SpliceGrapher/formats/parsers/splice_graph.py tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_splice_graph.py
```

```bash
uv run mypy SpliceGrapher/formats/parsers/__init__.py SpliceGrapher/formats/parsers/splice_graph.py tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_splice_graph.py
```

**Step 2: Run the broader parser/gene-model slice**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model.py tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_splice_graph.py tests/test_integration_simple.py
```

Expected: pass.

---

### Task 4: Run full verification and commit

**Files:**
- Verify whole worktree

**Step 1: Run the full SGN suite**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider
```

**Step 2: Run the build**

Run:

```bash
uv build
```

**Step 3: Inspect the final diff**

Run:

```bash
git status --short
git diff --stat main...HEAD
```

**Step 4: Commit**

Use a scoped message such as:

```bash
git commit -m "refactor(parsers): tighten facade and modernize splice graph parser"
```
