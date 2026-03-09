# Replace `formats.models` Flat Module With Package Directory Design

## Context

After `#175`, `SpliceGrapher/formats/models.py` is no longer a domain monolith. It is a 119-line compatibility aggregator that re-exports symbols from two real implementation modules:

- `SpliceGrapher/formats/model_domain.py`
- `SpliceGrapher/formats/model_index.py`

That split paid down the internal complexity, but it left the package layout in an intermediate state:

- the parent `formats/` directory still carries a flat `models.py` shell
- the internal implementation files still have transitional `model_*` prefixes
- one module-layout test still imports the temporary top-level `model_domain` and `model_index` modules directly

The runtime import boundary we actually want to preserve is `SpliceGrapher.formats.models`, not the file layout that produced it.

## Goal

Replace the flat `SpliceGrapher/formats/models.py` file with a real package directory:

- `SpliceGrapher/formats/models/__init__.py`
- `SpliceGrapher/formats/models/domain.py`
- `SpliceGrapher/formats/models/index.py`

The public import surface must remain stable for runtime callers using `SpliceGrapher.formats.models`, while the temporary top-level modules `model_domain.py` and `model_index.py` are deleted.

## Options Considered

### Option 1: Keep `formats/models.py` as a flat aggregator

This has the smallest diff, but it preserves the transitional layout debt indefinitely. The package boundary remains conceptually wrong because the real implementation is already multi-file.

### Option 2: Delete `formats.models` entirely and force callers onto `model_domain` / `model_index`

This is a hard cut, but it is the wrong hard cut. Runtime callers conceptually depend on a `models` namespace, and deleting that name would create needless import churn without improving the actual internal layout.

### Option 3: Replace `formats/models.py` with a real `formats/models/` package

This keeps the stable namespace that callers already use while deleting the transitional top-level module layout. It also allows us to drop the `model_` filename prefixes and make the directory structure match the responsibility boundary.

## Decision

Use Option 3.

We will keep `SpliceGrapher.formats.models` as the public namespace, but it will become a package directory instead of a flat file. The implementation files will move under that directory and be renamed to:

- `SpliceGrapher/formats/models/domain.py`
- `SpliceGrapher/formats/models/index.py`

`SpliceGrapher/formats/models/__init__.py` will become the new aggregation surface with the existing `__all__`-driven re-export contract.

The transitional top-level files:

- `SpliceGrapher/formats/model_domain.py`
- `SpliceGrapher/formats/model_index.py`
- `SpliceGrapher/formats/models.py`

will all be deleted.

## Import Boundary

### Preserved

These imports must continue to work after the cut:

- `from SpliceGrapher.formats import models`
- `import SpliceGrapher.formats.models as model_domain`
- `from SpliceGrapher.formats.models import Gene`
- all current runtime uses in:
  - `SpliceGrapher/formats/gene_model.py`
  - `SpliceGrapher/formats/serializers.py`
  - `SpliceGrapher/formats/parsers/gene_model_gff.py`
  - `SpliceGrapher/formats/parsers/gene_model_gff_context.py`
  - `SpliceGrapher/formats/parsers/gene_model_gff_records.py`
  - `SpliceGrapher/formats/parsers/gene_model_gff_handlers.py`

### Deleted

These transitional imports will be removed and any remaining callers updated:

- `from SpliceGrapher.formats import model_domain`
- `from SpliceGrapher.formats import model_index`
- `import SpliceGrapher.formats.model_domain`
- `import SpliceGrapher.formats.model_index`

Current evidence shows those top-level names are no longer used by runtime code and only remain in one test module.

## Internal Structure

### `SpliceGrapher/formats/models/__init__.py`

Responsibilities:

- re-export the public domain and index symbols
- own the public `__all__`
- preserve the stable `SpliceGrapher.formats.models` namespace

### `SpliceGrapher/formats/models/domain.py`

Responsibilities:

- own gene/transcript/feature domain entities
- own domain-level constants and helper functions currently in the former `model_domain.py`
- import index helpers via relative package imports where required

### `SpliceGrapher/formats/models/index.py`

Responsibilities:

- own chromosome/index helpers and interval query functions currently in the former `model_index.py`
- use `TYPE_CHECKING` and postponed annotations to avoid runtime cycles back into `domain.py`

## Testing Strategy

We do not need a new behavioral feature. We need a filesystem/layout hard cut with import stability.

Tests should prove two things:

1. runtime import stability for `SpliceGrapher.formats.models`
2. deletion of the temporary top-level `model_domain` / `model_index` modules

The existing module-layout test should be rewritten to:

- import `SpliceGrapher.formats.models` as the stable public boundary
- import the new package internals via:
  - `import SpliceGrapher.formats.models.domain as model_domain`
  - `import SpliceGrapher.formats.models.index as model_index`
- assert the aggregator still re-exports the same symbols from the new package internals

Gene-model parser, serializer, and integration tests should remain the regression net for runtime behavior.

## Risks

### Risk: circular imports reappear after moving files under one package

Mitigation:

- keep `from __future__ import annotations`
- preserve or tighten `TYPE_CHECKING` imports in the index layer
- prefer relative imports inside the package so the dependency graph stays explicit

### Risk: tests or docs still refer to the deleted top-level modules

Mitigation:

- search the repo for `model_domain` and `model_index`
- update the module-layout test and any plan/invariant text that asserts the older transitional layout

### Risk: package discovery/import semantics differ from the flat file

Mitigation:

- keep `__init__.py` explicit
- rerun mypy, targeted gene-model tests, full SGN `pytest`, and `uv build`

## Out of Scope

- changing the public `SpliceGrapher.formats.models` namespace
- further domain refactors inside `domain.py` or `index.py`
- backend parser changes such as `polars-bio` or `biobear` adoption
- additional gene-model API cleanup outside the package-layout cut
