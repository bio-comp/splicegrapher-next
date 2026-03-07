# Alignment IO Contract Cleanup Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Harden `SpliceGrapher/formats/alignment_io.py` tuple-shape guarantees and internal structure without changing its public ABI.

**Architecture:** Keep public camelCase entrypoints and tuple return types stable. Add failing tests for the `alignments=True` contract and depths-file fallback, introduce a minimal protocol for `pysamReadDepths`, and reorganize internal helpers into source-opening, pysam-collection, fallback-bridge, and public-wrapper blocks.

**Tech Stack:** Python 3.13, pytest, mypy, Ruff, pysam, numpy.

---

### Task 1: Lock The Tuple-Shape Contract In Tests

**Files:**
- Modify: `tests/test_splicegrapher_alignment_io.py`
- Modify: `tests/test_alignment_io_process_utils_boundary.py`
- Modify: `tests/test_alignment_io_parity.py`

**Step 1: Write the failing BAM/SAM alignment-map test**

Add a focused test that exercises the direct pysam path:

```python
def test_get_sam_read_data_alignments_true_returns_alignment_map(tmp_path: Path) -> None:
    fixture = build_alignment_fixture(tmp_path, repeat_scale=8)

    depths, junctions, alignments = getSamReadData(str(fixture.bam), alignments=True)

    assert fixture.chrom in depths
    assert fixture.chrom in junctions
    assert fixture.chrom in alignments
    assert alignments[fixture.chrom]
```

**Step 2: Run the single test to verify current behavior**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_splicegrapher_alignment_io.py -k alignments_true_returns_alignment_map
```

Expected: PASS or near-pass once the test is wired correctly.

**Step 3: Write the failing depths-fallback tuple-shape test**

Add a monkeypatched fallback test in `tests/test_alignment_io_process_utils_boundary.py`:

```python
def test_get_sam_read_data_depths_fallback_preserves_three_tuple_when_alignments_requested(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    sample_path = tmp_path / "sample.depths"
    sample_path.write_text("placeholder\n", encoding="utf-8")

    fake_depths = {"chr1": [0, 1, 2]}
    fake_junctions = {"chr1": []}

    monkeypatch.setattr(alignment_io, "_is_depths_source", lambda source: True)
    monkeypatch.setattr(
        alignment_io,
        "read_depths",
        lambda *args, **kwargs: (fake_depths, fake_junctions),
    )

    depths, junctions, alignments = getSamReadData(str(sample_path), alignments=True)

    assert list(depths["chr1"]) == [0, 1, 2]
    assert junctions == fake_junctions
    assert alignments == {}
```

**Step 4: Run the focused fallback test and confirm it fails**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_process_utils_boundary.py -k preserves_three_tuple
```

Expected: FAIL with unpacking mismatch until the fallback is fixed.

**Step 5: Commit the red tests checkpoint**

```bash
git add tests/test_splicegrapher_alignment_io.py tests/test_alignment_io_process_utils_boundary.py tests/test_alignment_io_parity.py
git commit -m "test: pin alignment io tuple-shape contract"
```

### Task 2: Tighten Touched Types Without Changing Public ABI

**Files:**
- Modify: `SpliceGrapher/formats/alignment_io.py`
- Test: `tests/test_alignment_io_process_utils_boundary.py`

**Step 1: Add a structural protocol for read-depth gene bounds**

Add a minimal protocol near the type aliases:

```python
class GeneBounds(Protocol):
    id: str
    strand: str
    minpos: int
    maxpos: int
```

**Step 2: Update `pysamReadDepths` to use the protocol**

```python
def pysamReadDepths(
    bamFile: pysam.AlignmentFile,
    chromosome: str,
    gene: GeneBounds,
    *,
    margin: int = 0,
    verbose: bool = False,
) -> tuple[int, numpy.ndarray]:
    ...
```

**Step 3: Introduce clearer source aliases where they improve readability**

Prefer a small alias set such as:

```python
SamTextSource = io.IOBase | list[str] | tuple[str, ...] | set[str] | list[bytes] | tuple[bytes, ...] | set[bytes]
AlignmentSource = str | PathLike[str] | SamTextSource
```

Apply these only where they reduce ambiguity; do not churn the whole file for cosmetic typing.

**Step 4: Run MyPy on the touched module/tests**

Run:
```bash
uv run mypy SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py
```

Expected: PASS.

**Step 5: Commit the typing cleanup**

```bash
git add SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py
git commit -m "refactor: tighten alignment io typing"
```

### Task 3: Fix Depths Fallback Contract Drift

**Files:**
- Modify: `SpliceGrapher/formats/alignment_io.py`
- Test: `tests/test_alignment_io_process_utils_boundary.py`
- Test: `tests/test_splicegrapher_alignment_io.py`

**Step 1: Extract a depths-fallback bridge helper**

Add a helper that owns the `read_depths(...)` bridge and tuple shaping:

```python
def _collect_depths_source_data(
    source: DepthSource,
    *,
    alignments: bool,
    junctions: bool,
    maxpos: int,
    minanchor: int,
    minjct: int,
    verbose: bool,
) -> CollectResult | CollectResultWithAlignments:
    depths, jcts = read_depths(...)
    normalized_depths = _depth_map_to_arrays(depths)
    if alignments:
        return normalized_depths, jcts, {}
    return normalized_depths, jcts
```

**Step 2: Route `getSamReadData` through the helper**

Replace the inline fallback block with the helper call so the tuple shaping lives in one place.

**Step 3: Re-run the focused fallback tests**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_process_utils_boundary.py -k preserves_three_tuple
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_splicegrapher_alignment_io.py -k alignments_true_returns_alignment_map
```

Expected: PASS.

**Step 4: Run the full alignment slice**

Run:
```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS.

**Step 5: Commit the contract fix**

```bash
git add SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
git commit -m "fix: preserve alignment io tuple contract"
```

### Task 4: Reorganize Internal Helper Blocks For Future Splits

**Files:**
- Modify: `SpliceGrapher/formats/alignment_io.py`

**Step 1: Reorder helpers into four concern blocks**

Organize the file as:
1. type aliases and constants
2. source normalization/opening helpers
3. pysam collection helpers
4. depths fallback bridge
5. public ABI wrappers

**Step 2: Add minimal section comments only where they reduce ambiguity**

Example:

```python
# Depths-file fallback bridge for legacy public wrappers.
```

Avoid decorative comments.

**Step 3: Re-run Ruff format and Ruff check**

Run:
```bash
uv run ruff format SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
uv run ruff check SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS.

**Step 4: Re-run MyPy and the alignment pytest slice**

Run:
```bash
uv run mypy SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS.

**Step 5: Commit the internal reorganization**

```bash
git add SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
git commit -m "refactor: compartmentalize alignment io internals"
```

### Task 5: Final Branch Verification And PR Prep

**Files:**
- Modify: `docs/plans/2026-03-07-alignment-io-contract-cleanup-design.md`
- Modify: `docs/plans/2026-03-07-alignment-io-contract-cleanup-implementation-plan.md`

**Step 1: Record any design drift if implementation differs from plan**

Update the design/plan docs only if the landed code diverges materially.

**Step 2: Run final touched verification**

Run:
```bash
uv run ruff check SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
uv run ruff format --check SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
uv run mypy SpliceGrapher/formats/alignment_io.py tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_alignment_io_parity.py tests/test_alignment_io_process_utils_boundary.py tests/test_splicegrapher_alignment_io.py tests/test_shortread_compat.py
```

Expected: PASS.

**Step 3: Commit any final doc adjustments**

```bash
git add docs/plans/2026-03-07-alignment-io-contract-cleanup-design.md docs/plans/2026-03-07-alignment-io-contract-cleanup-implementation-plan.md
git commit -m "docs: finalize alignment io cleanup plan"
```
