# PROVENANCE.md

This document tracks code/data/document provenance for `splicegrapher-next`
during extraction from iDiffIR's vendored `SpliceGrapher` tree and any future
upstream/downstream integration work.

## Repository Identity

- Repository: `bio-comp/splicegrapher-next`
- Purpose: generally usable Python 3 modernization/continuation of
  SpliceGrapher, with compatibility priorities for iDiffIR and TAPIS
- Import namespace compatibility target: `SpliceGrapher`

## Current Source Inventory (Phase: Active Extraction)

Current tracked extraction baseline includes:

- Identity/package artifacts: `SpliceGrapher/__init__.py`
- Core/shared baseline: `SpliceGrapher/SpliceGraph.py`,
  `SpliceGrapher/shared/*`
- Formats slice A baseline: `SpliceGrapher/formats/*` (current extracted files)
- Residual batch-0 identity completion: `SpliceGrapher/SpliceGrapher.cfg`

## Upstream and Derived Sources

Planned/known source relationships:

1. Historical upstream SpliceGrapher project
   - URL: <https://splicegrapher.sourceforge.net/>
   - Publication reference:
     Rogers MF, Thomas J, Reddy AS, Ben-Hur A. Genome Biol. 2012;13(1):R4.
2. iDiffIR vendored/refactored SpliceGrapher tree (primary extraction source)
   - Local source path used for extraction work:
     `iDiffIR/iDiffIR/SpliceGrapher`
   - Repository context: <https://github.com/bio-comp/iDiffIR>
3. TAPIS compatibility consumer context
   - Repository context: <https://github.com/bio-comp/tapis>

## Licensing Notes (Non-Legal Summary)

- Top-level license currently present in this repository: `LICENSE`.
- Licensing/provenance source of truth is centralized in project-level docs:
  `LICENSE`, `NOTICE`, `PROVENANCE.md`, and `LICENSES/`.
- Vendored iDiffIR SpliceGrapher tree includes mixed licensing provenance;
  required notices must remain discoverable through project inventory docs.
- Do not treat this file as legal advice; it is an engineering traceability
  record.

## Attribution and Notice Handling Rules

When copying or normalizing files in this repository:

- Preserve attribution and provenance through project-level inventories.
- Remove redundant per-file legacy license header blocks in tracked cleanup
  issues when legally safe.
- Preserve file-level notices when explicitly required by source licensing
  terms, and document those exceptions here.
- Record source path, commit/ref, and destination path in this file.
- Document any intentional modifications from the copied baseline.

## Extraction Ledger Template

Use one entry per extraction change-set:

```text
Date:
Issue/PR:
Source repository/path:
Source commit/ref:
Destination paths:
License/notice files added or updated:
Behavioral modifications from source baseline:
Validation run:
```

## Initial Follow-Up Items

- Add explicit source commit hashes for each copied batch.
- Add file-level license inventory rows in `LICENSES/README.md`.

## Extraction Ledger Entries

### 2026-02-22 - Core/shared baseline extraction (issue #5)

Date:

- 2026-02-22

Issue/PR:

- Issue: #5
- PR: (to be filled when opened)

Source repository/path:

- Repository: `bio-comp/iDiffIR`
- Path: `iDiffIR/SpliceGrapher`

Source commit/ref:

- `e69b8ada59860d270f6cecb2468d5993a995fac8`

Destination paths:

- `SpliceGrapher/__init__.py`
- `SpliceGrapher/SpliceGraph.py`
- `SpliceGrapher/shared/__init__.py`
- `SpliceGrapher/shared/utils.py`
- `SpliceGrapher/shared/streams.py`
- `tests/test_core_shared_import_smoke.py`

License/notice files added or updated:

- None in this extraction slice.

Behavioral modifications from source baseline:

- Updated internal import namespace from `iDiffIR.SpliceGrapher.*` to
  `SpliceGrapher.*` for extracted modules.
- Added import smoke tests for extracted core/shared modules.

Validation run:

- `uv run pytest tests/test_core_shared_import_smoke.py -q`
- `uv run python -c "import SpliceGrapher, SpliceGrapher.SpliceGraph, SpliceGrapher.shared.utils"`
- `uv build`

### 2026-02-24 - Residual identity cfg extraction (issue #6)

Date:

- 2026-02-24

Issue/PR:

- Issue: #6
- PR: (to be filled when opened)

Source repository/path:

- Repository: `bio-comp/iDiffIR`
- Path: `iDiffIR/SpliceGrapher/SpliceGrapher.cfg`

Source commit/ref:

- `e69b8ada59860d270f6cecb2468d5993a995fac8`

Destination paths:

- `SpliceGrapher/SpliceGrapher.cfg`

License/notice files added or updated:

- None in this extraction slice.

Behavioral modifications from source baseline:

- None. File copied as-is to complete batch-0 identity compatibility artifact.

Validation run:

- `uv run python -c "import SpliceGrapher"`
- `uv build`

### 2026-02-27 - License-header normalization tranche A (issue #34)

Date:

- 2026-02-27

Issue/PR:

- Issue: #34
- PR: (to be filled when opened)

Source repository/path:

- Repository: `bio-comp/splicegrapher-next`
- Paths:
  - `SpliceGrapher/SpliceGraph.py`
  - `SpliceGrapher/shared/streams.py`
  - `SpliceGrapher/shared/ShortRead.py`
  - `SpliceGrapher/formats/fasta.py`
  - `SpliceGrapher/formats/GeneModel.py`
  - `SpliceGrapher/formats/loader.py`
  - `SpliceGrapher/formats/alignment_io.py`

Source commit/ref:

- `c1ea06b1a8f7fdb5f0bff8984180ef6672a8ea30`

Destination paths:

- same as source paths above (in-place normalization)

License/notice files added or updated:

- `NOTICE`
- `PROVENANCE.md`

Behavioral modifications from source baseline:

- None. Removed redundant top-of-file legacy license header comment blocks.
- Standardized repository policy to project-level licensing/provenance
  inventory docs with exception handling for required file-level notices.

Validation run:

- `uv run ruff check . --fix`
- `uv run ruff format .`
- `uv run pytest`
- `uv build`

### 2026-03-01 - Remove shared streams module (issue #58)

Date:

- 2026-03-01

Issue/PR:

- Issue: #58
- PR: (to be filled when opened)

Source repository/path:

- Repository: `bio-comp/splicegrapher-next`
- Paths:
  - `SpliceGrapher/shared/process_utils.py`
  - `SpliceGrapher/shared/streams.py`
  - `tests/test_process_utils.py`
  - `tests/test_core_shared_import_smoke.py`

Source commit/ref:

- `7158bf0ec8f5e8ee761ad337839c2f7b3958f8c9`

Destination paths:

- `SpliceGrapher/shared/process_utils.py`
- `SpliceGrapher/shared/streams.py` (removed)
- `tests/test_process_utils.py`
- `tests/test_core_shared_import_smoke.py`

License/notice files added or updated:

- `PROVENANCE.md`

Behavioral modifications from source baseline:

- Replaced process-global fd redirection in `runCommand` with subprocess-local
  redirection using `subprocess.DEVNULL` defaults.
- Removed `SpliceGrapher/shared/streams.py` and all code/test imports of
  `SpliceGrapher.shared.streams`.

Validation run:

- `uv run pytest -q tests/test_process_utils.py tests/test_core_shared_import_smoke.py`
- `uv run ruff check .`
- `uv run ruff format --check .`
- `uv run pytest -q`
- `uv build`
