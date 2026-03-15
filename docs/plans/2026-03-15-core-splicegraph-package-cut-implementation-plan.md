# Core SpliceGraph Package Cut Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the flat `SpliceGrapher/core/splice_graph.py` module with a concrete package split across `constants.py`, `node.py`, and `graph.py`, then update SGN imports and invariant guards to the new layout.

**Architecture:** This is a structural hard cut, not a semantic rewrite. The existing splice-graph code is redistributed into focused modules inside `SpliceGrapher/core/splice_graph/`, the old flat file is removed, and SGN callers are updated directly to concrete modules instead of a broad compatibility facade.

**Tech Stack:** Python 3.11+, `networkx`, `pytest`, `ruff`, `mypy`, Hatchling/`uv`

---

### Task 1: Pin the New Package Boundary in Tests

**Files:**
- Modify: `tests/test_splice_graph.py`
- Modify: `tests/test_clean_invariant_guard.py`
- Test: `tests/test_splice_graph.py`
- Test: `tests/test_clean_invariant_guard.py`

**Step 1: Write the failing test**

Add assertions that expect the package layout instead of the flat file:
- imports resolve from `SpliceGrapher.core.splice_graph.graph`
- imports resolve from `SpliceGrapher.core.splice_graph.node`
- imports resolve from `SpliceGrapher.core.splice_graph.constants`
- clean-invariant guard now protects the package module files, not `SpliceGrapher/core/splice_graph.py`

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_splice_graph.py tests/test_clean_invariant_guard.py`

Expected: FAIL because the package modules do not exist yet and invariant tests still reference the flat file.

**Step 3: Write minimal implementation**

Do not touch production code yet. Only land the failing expectations in the tests.

**Step 4: Run test to verify it fails cleanly**

Run the same command again and confirm the failure is specifically missing-module / stale-path related.

**Step 5: Commit**

```bash
git add tests/test_splice_graph.py tests/test_clean_invariant_guard.py
git commit -m "test: pin core splicegraph package boundary"
```

### Task 2: Create the Package Modules and Remove the Flat File

**Files:**
- Delete: `SpliceGrapher/core/splice_graph.py`
- Create: `SpliceGrapher/core/splice_graph/__init__.py`
- Create: `SpliceGrapher/core/splice_graph/constants.py`
- Create: `SpliceGrapher/core/splice_graph/node.py`
- Create: `SpliceGrapher/core/splice_graph/graph.py`
- Test: `tests/test_splice_graph.py`

**Step 1: Write the failing test**

Use the new package-boundary test from Task 1 as the red test.

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_splice_graph.py`

Expected: FAIL because `graph.py`, `node.py`, and `constants.py` are not yet implemented.

**Step 3: Write minimal implementation**

Move code out of the flat file into focused modules:
- constants and type alias into `constants.py`
- `NullNode` and `SpliceGraphNode` into `node.py`
- `SpliceGraph` into `graph.py`
- keep `__init__.py` minimal

**Step 4: Run test to verify it passes**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_splice_graph.py`

Expected: PASS.

**Step 5: Commit**

```bash
git add SpliceGrapher/core/splice_graph tests/test_splice_graph.py
git commit -m "refactor(core): split splicegraph core package"
```

### Task 3: Rewrite SGN Imports and Invariant Guards

**Files:**
- Modify: `SpliceGrapher/core/graph_math.py`
- Modify: `SpliceGrapher/core/splicing_events.py`
- Modify: `SpliceGrapher/formats/parsers/splice_graph.py`
- Modify: `SpliceGrapher/formats/writers/splice_graph.py`
- Modify: `scripts/ci/check_clean_invariant.py`
- Modify: `tests/test_graph_math.py`
- Modify: `tests/test_integration_simple.py`
- Modify: `tests/test_splicing_events.py`
- Modify: `tests/test_clean_invariant_guard.py`

**Step 1: Write the failing test**

Keep the focused slice red by updating only one import at a time until all stale imports/path assumptions are visible.

**Step 2: Run test to verify it fails**

Run: `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_graph_math.py tests/test_splicing_events.py tests/test_integration_simple.py tests/test_clean_invariant_guard.py`

Expected: FAIL on stale imports from `SpliceGrapher.core.splice_graph` or stale protected flat-file paths.

**Step 3: Write minimal implementation**

Update SGN imports directly to concrete modules:
- `graph.py` for `SpliceGraph`
- `node.py` for `SpliceGraphNode`
- `constants.py` for shared constants

Update invariant guard logic to protect the new package module files instead of the deleted flat file.

**Step 4: Run test to verify it passes**

Run the same command again.

Expected: PASS.

**Step 5: Commit**

```bash
git add SpliceGrapher/core/graph_math.py SpliceGrapher/core/splicing_events.py SpliceGrapher/formats/parsers/splice_graph.py SpliceGrapher/formats/writers/splice_graph.py scripts/ci/check_clean_invariant.py tests/test_graph_math.py tests/test_integration_simple.py tests/test_splicing_events.py tests/test_clean_invariant_guard.py
git commit -m "refactor: rewrite splicegraph imports to concrete modules"
```

### Task 4: Update Current-State Docs and Run Full Verification

**Files:**
- Modify: `PROVENANCE.md`
- Modify: `tests/test_clean_invariant_guard.py`
- Test: repo-wide verification commands

**Step 1: Write the failing test**

No new behavior test is needed here; this task is for current-state doc accuracy and full verification.

**Step 2: Run verification before doc touch**

Run:
- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run mypy .`
- `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
- `uv build`

Expected: PASS before final docs commit, or reveal any missed import/layout fallout.

**Step 3: Write minimal implementation**

Update `PROVENANCE.md` current-source inventory language so it references the package layout instead of the deleted flat file.

**Step 4: Run verification to verify it still passes**

Run the full command set again if `ruff format` or doc edits touched tracked files in ways that need re-checking.

**Step 5: Commit**

```bash
git add PROVENANCE.md
git commit -m "docs: update splicegraph package provenance"
```
