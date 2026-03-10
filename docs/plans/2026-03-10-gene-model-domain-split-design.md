# Gene-Model Domain Split Design

Issue: #191  
Branch: `sg-next/191-model-domain-split`

## Problem

`SpliceGrapher/formats/models/domain.py` is now the largest remaining SGN runtime monolith. It mixes:
- shared constants and factories
- the `Locus` value object
- base feature classes
- transcript assembly logic
- gene assembly logic

That shape is still import-stable through `SpliceGrapher.formats.models`, but it is not internally maintainable.

## Goal

Split `SpliceGrapher/formats/models/domain.py` into smaller internal modules while preserving the `SpliceGrapher.formats.models` package boundary.

## Non-Goals

- no caller-facing import break outside `SpliceGrapher.formats.models`
- no parser behavior changes
- no `annotation_io` or parser-handler refactors in this slice

## Target Shape

`SpliceGrapher/formats/models/`
- `__init__.py`
- `constants.py`
- `locus.py`
- `features.py`
- `transcript.py`
- `gene.py`
- `domain.py` removed
- `index.py` kept in place

## Constraints

- keep `from SpliceGrapher.formats.models import Gene` stable
- avoid circular imports with `from __future__ import annotations` and `TYPE_CHECKING`
- keep tests focused on layout plus behavior parity

## Verification

- targeted gene-model test slice
- full `pytest`
- `ruff`, `mypy`, `uv build`
