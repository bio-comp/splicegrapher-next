# Gene Model Record Handlers Hard Clean Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Hard-clean `gene_model_gff_record_handlers.py` by extracting internal dispatch helpers while preserving parser behavior.

**Architecture:** Keep handler entrypoints and `RECORD_HANDLERS` unchanged, but move repeated transcript-parent resolution and transcript-link registration into private helper seams used by exon, mRNA, and transcript-region handlers.

**Tech Stack:** Python 3.11+, pytest, ruff, mypy, uv.

---

### Task 1: Add red boundary tests for helper seams

**Files:**
- Create: `tests/test_gene_model_gff_record_handlers_boundary.py`

1. Add tests requiring private helpers for transcript-link registration and region-target resolution.
2. Run red test command:
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_gff_record_handlers_boundary.py`
3. Verify expected failure is missing helper attributes.

### Task 2: Extract helper seams in handlers module

**Files:**
- Modify: `SpliceGrapher/formats/parsers/gene_model_gff_record_handlers.py`

1. Add `_parent_gene_for_transcript(...)` for centralized transcript->gene lookup.
2. Add `_register_transcript_links(...)` for model map registration and pending drain.
3. Add `_resolve_region_target(...)` for CDS/UTR parent resolution with pending queue behavior.
4. Rewire `handle_exon_record`, `handle_mrna_record`, and `handle_transcript_region_record` to use helpers.

### Task 3: Verify parity and run full gates

**Files:**
- Verify: parser and gene-model tests

1. Run focused parser/gene-model suite:
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_gff_record_handlers_boundary.py tests/test_gene_model_gff_parser.py tests/test_gene_model_gff_module_layout.py tests/test_gene_model.py`
2. Run full gate chain:
   - `uv run ruff check . --fix`
   - `uv run ruff format .`
   - `uv run mypy .`
   - `PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider`
   - `uv build`
3. Commit only tracked code/test/docs changes and open PR linked to `#210`.
