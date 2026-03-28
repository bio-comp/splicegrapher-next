# FASTA Backend Pyfastx Evaluation Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Produce an evidence-backed pyfastx backend recommendation while keeping SGN runtime FASTA behavior unchanged.

**Architecture:** Add an isolated benchmark probe and decision artifacts only; do not touch runtime `SpliceGrapher.formats.fasta` behavior in this issue.

**Tech Stack:** Python 3.11+, pyfaidx, optional pyfastx (`uv run --with pyfastx`), pytest, ruff, mypy.

---

### Task 1: Add isolated backend probe

**Files:**
- Create: `benchmarks/fasta_backend_probe.py`

1. Implement SGN/pyfaidx/pyfastx parity + performance probe.
2. Include contract checks and indexing overhead capture.
3. Keep pyfastx optional and behavior-neutral for runtime code.

### Task 2: Add probe regression test (no pyfastx required)

**Files:**
- Create: `tests/test_fasta_backend_probe.py`

1. Validate probe output shape and baseline parity using pyfaidx + SGN path.
2. Force pyfastx-absent path via monkeypatch and verify recommendation state.

### Task 3: Generate artifact + publish decision

**Files:**
- Create: `docs/testing/2026-03-25-fasta-backend-probe.json`
- Create: `docs/adr/2026-03-25-fasta-backend-pyfastx-evaluation.md`

1. Run probe using:
   - `uv run --with pyfastx python benchmarks/fasta_backend_probe.py --output docs/testing/2026-03-25-fasta-backend-probe.json`
2. Write ADR with recommendation and downstream compatibility implications.

### Task 4: Verify and prepare PR

1. Run focused tests:
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_fasta_io.py tests/test_fasta_module_layout.py tests/test_fasta_backend_probe.py`
2. Run full quality gates:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
3. Open PR linked to `#53` with recommendation summary.
