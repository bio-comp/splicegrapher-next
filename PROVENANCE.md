# PROVENANCE.md

This document tracks code/data/document provenance for `splicegrapher-next`
during extraction from iDiffIR's vendored `SpliceGrapher` tree and any future
upstream/downstream integration work.

## Repository Identity

- Repository: `bio-comp/splicegrapher-next`
- Purpose: shared Python 3 modernization/continuation bridge for SpliceGrapher
  used by iDiffIR and TAPIS
- Import namespace compatibility target: `SpliceGrapher`

## Current Source Inventory (Phase: Scaffold)

This repository currently contains bootstrap packaging and documentation assets.
No large SpliceGrapher source-module extraction is recorded yet in this file.

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
- Vendored iDiffIR SpliceGrapher tree includes additional file-level licensing
  context that must be preserved during extraction.
- Do not treat this file as legal advice; it is an engineering traceability
  record.

## Attribution and Header Preservation Rules

When copying files into this repository:

- Preserve original copyright headers.
- Preserve file-level notices and author attributions.
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

- Populate this ledger when issue #6 extraction begins.
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
- `uv run python -c "import SpliceGrapher, SpliceGrapher.SpliceGraph, SpliceGrapher.shared.utils, SpliceGrapher.shared.streams"`
- `uv build`
