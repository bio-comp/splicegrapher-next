# Alignment IO Hard-Break Rewrite Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the legacy camelCase alignment I/O API with a snake_case-only surface and repair SGN call sites/tests to the new boundary.

**Architecture:** Rewrite tests first to the new API, then rename/delete the public functions inside `SpliceGrapher/formats/alignment_io.py`, preserving the corrected tuple-shape behavior in `collect_alignment_data(...)`. Use failing imports/calls as the migration map across SGN.

**Tech Stack:** Python 3.13, pytest, mypy, Ruff, pysam, numpy.

---

### Task 1: Rewrite Alignment Tests To The New Public API

**Files:**
- Modify: `tests/test_alignment_io_parity.py`
- Modify: `tests/test_alignment_io_process_utils_boundary.py`
- Modify: `tests/test_splicegrapher_alignment_io.py`

**Step 1: Replace camelCase imports with snake_case imports**

Change imports to the new names:
- `collect_alignment_data`
- `read_alignment_depths`
- `read_alignment_junctions`
- `read_alignment_spans`
- `read_alignment_headers`

**Step 2: Rewrite test call sites to the new API**

Examples:
```python
sam_depths, sam_junctions = collect_alignment_data(str(fixture.sam))
depths = read_alignment_depths(str(fixture.bam), chromosomes=[fixture.chrom])
junctions = read_alignment_junctions(str(fixture.bam), chromosomes=[fixture.chrom], minjct=1)
```

**Step 3: Run the focused alignment tests and confirm import/call failures**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: FAIL because the new API does not exist yet.

**Step 4: Commit the red test rewrite**

```bash
git add tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
git commit -m "test: rewrite alignment io tests to snake_case API"
```

### Task 2: Replace The Public alignment_io API

**Files:**
- Modify: `SpliceGrapher/formats/alignment_io.py`

**Step 1: Rename public entrypoints**

Rename the public API to:
- `collect_alignment_data`
- `read_alignment_depths`
- `read_alignment_junctions`
- `read_alignment_spans`
- `read_alignment_headers`
- `read_alignment_chromosome_info`
- `read_alignment_sequences`
- `calculate_gene_depths`

**Step 2: Rename helper functions that should not remain public**

Convert:
- `isBamFile` -> `_is_bam`
- `isCramFile` -> `_is_cram`
- `makeChromosomeSet` -> `_make_chromosome_set`
- `pysamStrand` -> `_record_strand_for_read` or inline where appropriate

**Step 3: Delete camelCase wrappers entirely**

Do not leave aliases or passthrough wrappers behind.

**Step 4: Keep the tuple-shape fix in the new API**

`collect_alignment_data(..., include_alignments=True)` must still return a 3-tuple for both the pysam path and the depths-file fallback path.

**Step 5: Run the focused alignment slice**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: either PASS or fail only on remaining SGN call sites/imports.

**Step 6: Commit the API cut**

```bash
git add SpliceGrapher/formats/alignment_io.py
git commit -m "refactor: replace alignment io public API"
```

### Task 3: Sweep SGN Call Sites To The New API

**Files:**
- Modify any SGN callers found by grep in `SpliceGrapher/`, `tests/`, `scripts/`, and `docs/`

**Step 1: Grep for removed names**

Run:
```bash
rg -n "getSamReadData|getSamDepths|getSamJunctions|getSamAlignments|getSamHeaders|getSamHeaderInfo|getSamSequences|pysamReadDepths|pysamStrand|isBamFile|isCramFile|makeChromosomeSet" SpliceGrapher tests scripts docs
```

**Step 2: Update each SGN call site to the new names**

Examples:
```python
collect_alignment_data(...)
read_alignment_depths(...)
read_alignment_junctions(...)
read_alignment_spans(...)
read_alignment_headers(...)
read_alignment_chromosome_info(...)
read_alignment_sequences(...)
calculate_gene_depths(...)
```

**Step 3: Re-run the focused alignment slice after each group of call-site edits**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS once the SGN call sites are updated.

**Step 4: Commit the caller migration**

```bash
git add SpliceGrapher tests scripts docs
git commit -m "refactor: migrate SGN to new alignment io API"
```

### Task 4: Final Touched Verification

**Files:**
- Modify docs/plans only if implementation drifted materially from the design

**Step 1: Run final touched verification**

Run:
```bash
uv run ruff check SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
uv run ruff format --check SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
uv run mypy SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS.

**Step 2: Commit any final doc adjustments**

```bash
git add docs/plans/2026-03-07-alignment-io-contract-cleanup-design.md docs/plans/2026-03-07-alignment-io-contract-cleanup-implementation-plan.md
git commit -m "docs: finalize alignment io hard-break plan"
```
