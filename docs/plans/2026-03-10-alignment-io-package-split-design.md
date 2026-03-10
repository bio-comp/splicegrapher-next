# Alignment I/O Package Split Design

Issue: #195  
Branch: `sg-next/195-alignment-io-package-split`

## Problem

`SpliceGrapher/formats/alignment_io.py` is now the largest remaining SGN runtime monolith. It still mixes:
- source/path normalization
- SAM/BAM/CRAM opening
- depth/junction/span collection
- depths-file fallback bridging
- public alignment loading APIs

The API is already modernized, but the internal file is still too large and mixed to maintain cleanly.

## Goal

Split `alignment_io.py` into a smaller internal package while preserving the current snake_case public API.

## Non-Goals

- no public API break for the current snake_case alignment functions
- no `depth_io` or `junction` semantic changes in this slice
- no extraction/package moves outside `SpliceGrapher/formats`

## Target Shape

`SpliceGrapher/formats/alignment_io/`
- `__init__.py`
- `types.py`
- `sources.py`
- `collect.py`
- `depths.py`
- `api.py`

## Constraints

- preserve current snake_case alignment public API
- keep tests focused on module layout plus behavior parity
- avoid reintroducing deleted camelCase compatibility wrappers

## Verification

- focused alignment/depth/shortread test slice
- full `pytest`
- `ruff`, `mypy`, `uv build`
