# SGN Enum Harmonization Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Harmonize SGN finite domains with canonical enums and migrate core control-flow paths to fail-fast enum usage, starting with SpliceGraph and GeneModel.

**Architecture:** Introduce central `str`-backed enums in `SpliceGrapher/core/enums.py`, add explicit boundary coercion helpers, then migrate core modules in small tranche commits with test-first changes. Internal logic fails fast on invalid enum values; serialization boundaries emit string `.value` to preserve output contracts unless intentionally changed.

**Tech Stack:** Python 3.11+, `Enum` (`str, Enum`), `pytest`, `ruff`, existing SGN CI clean-invariant guard.

---

### Task 1: Add Enum Foundation Module

**Files:**
- Create: `SpliceGrapher/core/enums.py`
- Create: `tests/test_core_enums.py`

**Step 1: Write the failing test**

```python
from SpliceGrapher.core.enums import Strand


def test_strand_values_are_stable() -> None:
    assert Strand.PLUS.value == "+"
    assert Strand.MINUS.value == "-"
    assert Strand.UNKNOWN.value == "."
```

**Step 2: Run test to verify it fails**

Run: `uv run pytest -q tests/test_core_enums.py`
Expected: FAIL with import error (`SpliceGrapher.core.enums` missing)

**Step 3: Write minimal implementation**

Implement `SpliceGrapher/core/enums.py` with:
- `Strand(str, Enum)`
- `RecordType(str, Enum)` (initial subset used by SpliceGraph + GeneModel)
- `NodeDisposition(str, Enum)`
- `EdgeType(str, Enum)`
- `AttrKey(str, Enum)`

**Step 4: Run test to verify it passes**

Run: `uv run pytest -q tests/test_core_enums.py`
Expected: PASS

**Step 5: Commit**

```bash
git add SpliceGrapher/core/enums.py tests/test_core_enums.py
git commit -m "feat(enums): add core enum foundation"
```

### Task 2: Add Boundary Coercion Helpers (Fail-Fast)

**Files:**
- Create: `SpliceGrapher/core/enum_coercion.py`
- Modify: `tests/test_core_enums.py`

**Step 1: Write the failing test**

```python
import pytest

from SpliceGrapher.core.enum_coercion import coerce_strand
from SpliceGrapher.core.enums import Strand


def test_coerce_strand_accepts_symbol() -> None:
    assert coerce_strand("+") is Strand.PLUS


def test_coerce_strand_rejects_invalid() -> None:
    with pytest.raises(ValueError):
        coerce_strand("?")
```

**Step 2: Run test to verify it fails**

Run: `uv run pytest -q tests/test_core_enums.py`
Expected: FAIL with import error for `enum_coercion`

**Step 3: Write minimal implementation**

Implement coercion helpers:
- `coerce_strand(value: str | Strand) -> Strand`
- `coerce_record_type(...) -> RecordType`
- `coerce_edge_type(...) -> EdgeType`

All helpers must raise `ValueError` on unknown input.

**Step 4: Run test to verify it passes**

Run: `uv run pytest -q tests/test_core_enums.py`
Expected: PASS

**Step 5: Commit**

```bash
git add SpliceGrapher/core/enum_coercion.py tests/test_core_enums.py
git commit -m "feat(enums): add fail-fast coercion helpers"
```

### Task 3: Add SpliceGraph Enum Regression Guard

**Files:**
- Modify: `tests/test_splice_graph.py`

**Step 1: Write the failing test**

Add test asserting invalid strand usage fails fast in enum-coerced path:

```python
import pytest

from SpliceGrapher.SpliceGraph import SpliceGraph


def test_union_rejects_invalid_other_strand_after_enum_migration() -> None:
    left = SpliceGraph("left", "chr1", "+")
    right = SpliceGraph("right", "chr1", "?")

    with pytest.raises(ValueError):
        left.union(right)
```

**Step 2: Run test to verify it fails**

Run: `uv run pytest -q tests/test_splice_graph.py`
Expected: FAIL (currently does not raise `ValueError`)

**Step 3: Write minimal implementation**

No implementation yet; this task only establishes red test before migration in Task 4.

**Step 4: Run test to keep red state visible**

Run: `uv run pytest -q tests/test_splice_graph.py`
Expected: FAIL

**Step 5: Commit**

```bash
git add tests/test_splice_graph.py
git commit -m "test(splicegraph): add fail-fast strand regression before enum migration"
```

### Task 4: Migrate SpliceGraph Core Domains to Enums

**Files:**
- Modify: `SpliceGrapher/SpliceGraph.py`
- Modify: `tests/test_splice_graph.py`

**Step 1: Implement minimal enum migration**

In `SpliceGrapher/SpliceGraph.py`:
- Import `Strand`, `RecordType`, `EdgeType`, `AttrKey`, coercion helpers.
- Replace string-domain comparisons in control flow with enum comparisons.
- Coerce incoming strand/type values at parser/boundary entry points.
- Preserve GFF text output by serializing enum values with `.value`.

**Step 2: Run focused tests**

Run: `uv run pytest -q tests/test_splice_graph.py`
Expected: PASS (including red test from Task 3)

**Step 3: Run file-level lint**

Run: `uv run ruff check SpliceGrapher/SpliceGraph.py tests/test_splice_graph.py`
Expected: PASS

**Step 4: Commit**

```bash
git add SpliceGrapher/SpliceGraph.py tests/test_splice_graph.py
git commit -m "refactor(splicegraph): migrate control-flow domains to enums"
```

### Task 5: Add GeneModel Enum Regression Tests

**Files:**
- Modify: `tests/test_gene_model.py`

**Step 1: Write the failing test**

Add a coercion-focused test that exercises gene/record type handling via enum-backed path.

```python
import pytest

from SpliceGrapher.core.enum_coercion import coerce_record_type


def test_coerce_record_type_rejects_unknown() -> None:
    with pytest.raises(ValueError):
        coerce_record_type("totally_unknown_type")
```

**Step 2: Run test to verify it fails**

Run: `uv run pytest -q tests/test_gene_model.py`
Expected: FAIL until GeneModel path integrates coercion in Task 6.

**Step 3: Keep test and red state**

No implementation yet.

**Step 4: Commit**

```bash
git add tests/test_gene_model.py
git commit -m "test(genemodel): add enum fail-fast regression coverage"
```

### Task 6: Migrate GeneModel Record Domains to Enums

**Files:**
- Modify: `SpliceGrapher/formats/GeneModel.py`
- Modify: `tests/test_gene_model.py`

**Step 1: Implement minimal enum migration**

In `SpliceGrapher/formats/GeneModel.py`:
- Replace record/strand string-domain control flow with enums/coercion.
- Keep serialization and external parsed values stable at boundaries.
- Use `.value` for text output contracts.

**Step 2: Run focused tests**

Run: `uv run pytest -q tests/test_gene_model.py`
Expected: PASS

**Step 3: Run lint for touched files**

Run: `uv run ruff check SpliceGrapher/formats/GeneModel.py tests/test_gene_model.py`
Expected: PASS

**Step 4: Commit**

```bash
git add SpliceGrapher/formats/GeneModel.py tests/test_gene_model.py
git commit -m "refactor(genemodel): migrate record and strand domains to enums"
```

### Task 7: Extend Clean-Invariant Guardrail for Magic Strings

**Files:**
- Modify: `scripts/ci/check_clean_invariant.py`
- Modify: `tests/test_clean_invariant_guard.py`

**Step 1: Write failing guard test**

Add unit test asserting new banned-pattern check flags magic-string control flow in protected modules.

**Step 2: Run test to verify it fails**

Run: `uv run pytest -q tests/test_clean_invariant_guard.py`
Expected: FAIL until guardrail logic is added.

**Step 3: Implement guardrail**

Add detection for banned magic-string patterns in:
- `SpliceGrapher/SpliceGraph.py`
- `SpliceGrapher/formats/GeneModel.py`

Allowlist boundary adapter modules only.

**Step 4: Run tests**

Run: `uv run pytest -q tests/test_clean_invariant_guard.py`
Expected: PASS

**Step 5: Commit**

```bash
git add scripts/ci/check_clean_invariant.py tests/test_clean_invariant_guard.py
git commit -m "ci: add enum migration guard against new magic-string control flow"
```

### Task 8: Update Tracking and Execute Full Verification

**Files:**
- Modify: `AGENTS.md` (local-only policy note)
- Modify: `LOCAL_TODO.md` (local tracker)

**Step 1: Update local policy/tracking**

- Add explicit enum harmonization directive and fail-fast posture.
- Add enum tranche checklist and active execution order.

**Step 2: Run full project verification**

Run:
- `uv run ruff check .`
- `uv run ruff format --check .`
- `uv run pytest -q`
- `uv run python scripts/ci/check_clean_invariant.py`

Expected: all PASS.

**Step 3: Commit local tracker/policy updates (if repository policy permits tracked changes)**

```bash
git add AGENTS.md LOCAL_TODO.md
git commit -m "docs(local): track enum harmonization policy and tranche checklist"
```

(If these files are untracked/local-only by policy, do not include this commit in PR scope.)

### Task 9: Open/Update Enum Epic Issue and Link Tranches

**Files:**
- Modify: issue tracker metadata (GitHub)

**Step 1: Create or update umbrella issue**

- Title: `Enum harmonization across SGN core domains (fail-fast migration)`
- Include tranche breakdown E0-E4, breakage policy, and verification gates.

**Step 2: Link tranche implementation issues**

- SpliceGraph tranche
- GeneModel/formats tranche
- shared/scripts tranche
- hardening tranche

**Step 3: Add PR linkage language**

- Ensure each PR references umbrella + tranche issue IDs and declares migration impact.

**Step 4: Commit**

No code commit; issue-tracker operation only.

