# Models Package Cut Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Replace the flat `SpliceGrapher/formats/models.py` compatibility file with a real `SpliceGrapher/formats/models/` package while preserving the `SpliceGrapher.formats.models` import surface.

**Architecture:** Move the current split internals into `SpliceGrapher/formats/models/{domain,index}.py`, convert `formats.models` into a package-level `__init__.py` aggregator, and delete the temporary top-level `model_domain.py` / `model_index.py` files. The cut is file-system structural, but tests must explicitly pin both the preserved public namespace and the deleted transitional import paths.

**Tech Stack:** Python 3.12, `uv`, `pytest`, `ruff`, `mypy`

---

### Task 1: Write the failing package-layout test

**Files:**
- Modify: `tests/test_gene_model_module_layout.py`

**Step 1: Rewrite the module-layout test to the target package shape**

Update the test so it imports:

```python
from SpliceGrapher.formats import models
import SpliceGrapher.formats.models.domain as model_domain
import SpliceGrapher.formats.models.index as model_index
```

Assert that:
- `models.Gene is model_domain.Gene`
- `models.Transcript is model_domain.Transcript`
- `models.BaseFeature is model_domain.BaseFeature`
- `models.ChromosomeGeneIndex is model_index.ChromosomeGeneIndex`
- `models.Chromosome is model_index.Chromosome`
- `models.feature_sort_key is model_index.feature_sort_key`

**Step 2: Run the test to verify it fails**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_module_layout.py
```

Expected: fail because `SpliceGrapher.formats.models` is still a flat module, not a package with `domain` / `index` submodules.

**Step 3: Commit the red test only after confirming failure**

Do not commit yet if the failure is not the expected package import failure.

---

### Task 2: Move the implementation files under a real package

**Files:**
- Create: `SpliceGrapher/formats/models/__init__.py`
- Create: `SpliceGrapher/formats/models/domain.py`
- Create: `SpliceGrapher/formats/models/index.py`
- Delete: `SpliceGrapher/formats/models.py`
- Delete: `SpliceGrapher/formats/model_domain.py`
- Delete: `SpliceGrapher/formats/model_index.py`

**Step 1: Create the package directory**

Create:

```text
SpliceGrapher/formats/models/
```

**Step 2: Move `model_domain.py` to `models/domain.py`**

Keep the implementation intact except for import rewrites.

**Step 3: Move `model_index.py` to `models/index.py`**

Keep the implementation intact except for import rewrites.

**Step 4: Convert `models.py` into `models/__init__.py`**

Move the aggregator content into `SpliceGrapher/formats/models/__init__.py` and rewrite imports to relative package imports:

```python
from .domain import ...
from .index import ...
```

Preserve `__all__`.

**Step 5: Delete the flat/transitional modules**

Delete:
- `SpliceGrapher/formats/models.py`
- `SpliceGrapher/formats/model_domain.py`
- `SpliceGrapher/formats/model_index.py`

**Step 6: Run the rewritten module-layout test**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_module_layout.py
```

Expected: pass.

**Step 7: Commit the package-layout cut**

```bash
git add SpliceGrapher/formats/models tests/test_gene_model_module_layout.py
git commit -m "refactor(formats): convert models module into package"
```

---

### Task 3: Repair internal imports and type-checking boundaries

**Files:**
- Modify: `SpliceGrapher/formats/gene_model.py`
- Modify: `SpliceGrapher/formats/serializers.py`
- Modify: `SpliceGrapher/formats/parsers/gene_model_gff.py`
- Modify: `SpliceGrapher/formats/parsers/gene_model_gff_context.py`
- Modify: `SpliceGrapher/formats/parsers/gene_model_gff_records.py`
- Modify: `SpliceGrapher/formats/parsers/gene_model_gff_handlers.py`

**Step 1: Update imports only where needed**

Preserve existing runtime calls that use:

```python
import SpliceGrapher.formats.models as model_domain
```

Only touch callers that reference deleted top-level modules or now-broken paths.

**Step 2: Repair internal package-relative imports**

Inside the new package files, replace old top-level imports such as:

```python
from SpliceGrapher.formats.model_domain import ...
from SpliceGrapher.formats.model_index import ...
```

with relative imports such as:

```python
from .domain import ...
from .index import ...
```

Use `TYPE_CHECKING` where runtime cycles would otherwise reappear.

**Step 3: Run targeted lint/type checks**

Run:

```bash
uv run ruff check SpliceGrapher/formats/models SpliceGrapher/formats/gene_model.py SpliceGrapher/formats/serializers.py SpliceGrapher/formats/parsers/gene_model_gff.py SpliceGrapher/formats/parsers/gene_model_gff_context.py SpliceGrapher/formats/parsers/gene_model_gff_records.py SpliceGrapher/formats/parsers/gene_model_gff_handlers.py tests/test_gene_model_module_layout.py
```

```bash
uv run mypy SpliceGrapher/formats/models SpliceGrapher/formats/gene_model.py SpliceGrapher/formats/serializers.py SpliceGrapher/formats/parsers/gene_model_gff.py SpliceGrapher/formats/parsers/gene_model_gff_context.py SpliceGrapher/formats/parsers/gene_model_gff_records.py SpliceGrapher/formats/parsers/gene_model_gff_handlers.py tests/test_gene_model_module_layout.py
```

Expected: both pass.

**Step 4: Commit import/typing repairs if they required separate edits**

Use a separate commit only if the package cut and import repair naturally split into reviewable changes.

---

### Task 4: Run the gene-model regression slice

**Files:**
- Test: `tests/test_gene_model_module_layout.py`
- Test: `tests/test_gene_model_gff_module_layout.py`
- Test: `tests/test_gene_model_gff_parser.py`
- Test: `tests/test_gene_model_serializers.py`
- Test: `tests/test_gene_model.py`
- Test: `tests/test_integration_simple.py`

**Step 1: Run the targeted gene-model regression slice**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider tests/test_gene_model_module_layout.py tests/test_gene_model_gff_module_layout.py tests/test_gene_model_gff_parser.py tests/test_gene_model_serializers.py tests/test_gene_model.py tests/test_integration_simple.py
```

Expected: pass.

**Step 2: If anything fails, fix import/layout fallout only**

Do not expand scope into additional gene-model refactors.

**Step 3: Commit regression fixes if needed**

Keep commit messages scoped to package import fallout, not unrelated cleanup.

---

### Task 5: Run full verification and open the PR

**Files:**
- Verify whole worktree only

**Step 1: Run the full SGN test suite**

Run:

```bash
PYTHONDONTWRITEBYTECODE=1 uv run pytest -q -p no:cacheprovider
```

Expected: pass.

**Step 2: Run the package build**

Run:

```bash
uv build
```

Expected: wheel and sdist build successfully.

**Step 3: Inspect the final diff**

Run:

```bash
git status --short
git diff --stat main...HEAD
```

Expected: only the package-layout cut and related test/import changes.

**Step 4: Commit any remaining tracked changes**

Use a final commit if needed.

**Step 5: Push and open the PR**

PR summary should state:
- `formats.models` is now a package directory, not a flat file
- top-level transitional `model_domain.py` / `model_index.py` are deleted
- runtime import surface remains `SpliceGrapher.formats.models`
