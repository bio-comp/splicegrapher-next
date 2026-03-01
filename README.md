# splicegrapher-next

`splicegrapher-next` is an unofficial Python 3 modernization/continuation of
SpliceGrapher for general splice-graph and RNA-seq analysis workflows, with
compatibility priorities for iDiffIR and TAPIS during active migration.

- Repository/distribution name: `splicegrapher-next`
- Python import namespace (compatibility): `SpliceGrapher`

## Why This Repo Exists

Historically, iDiffIR and TAPIS relied on divergent SpliceGrapher copies.
This repository reduces duplication, centralizes migration work, and provides
a generally usable Python 3 continuation while preserving behavior and imports.

## Phase 1 Scope

Goals:

- Extract and stabilize modernized SpliceGrapher components from iDiffIR.
- Preserve `SpliceGrapher` import compatibility where practical.
- Keep visualization (`view/`, `plot/`) as first-class functionality.
- Enable incremental external adoption by iDiffIR, TAPIS, and other users.

Non-goals:

- Full graph-engine rewrite.
- Big-bang API redesign.
- Immediate unification of iDiffIR/TAPIS project-specific logic.

## Quick Start

```bash
git clone git@github.com:bio-comp/splicegrapher-next.git
cd splicegrapher-next
```

### Development Install (uv-first)

```bash
uv sync --group dev
```

### Local Git Hook Gates (pre-commit)

Install local hooks:

```bash
uv run pre-commit install --hook-type pre-commit --hook-type pre-push
```

Manual hook runs:

```bash
uv run pre-commit run --all-files
uv run pre-commit run --hook-stage pre-push --all-files
```

Hook coverage:

- pre-commit stage: Ruff check/format + clean-invariant ratchet
- pre-push stage: pytest + `uv build`

### Editable pip Fallback

```bash
python -m pip install -e .
```

### Import Compatibility Check

```bash
uv run python -c "import SpliceGrapher; print(SpliceGrapher.__name__)"
```

## Configuration (Canonical TOML + Env Overrides)

`splicegrapher-next` now supports a canonical TOML configuration file validated
through `pydantic-settings`:

```toml
[splicegrapher]
gene_model = "/absolute/path/to/genes.gtf"
fasta_reference = "/absolute/path/to/reference.fa"
```

Use `splicegrapher.toml.example` as a starting point.

Configuration loading is TOML-only and intentionally does not support legacy
`.cfg` files.

Environment variables can override nested TOML values using:

- prefix: `SGN_`
- nested delimiter: `__`

Example override for orchestration layers (Nextflow/Docker env injection):

```bash
export SGN_SPLICEGRAPHER__GENE_MODEL=/mnt/data/genes.gtf
```

Detailed loader behavior and API examples are in `docs/configuration.md`.

## Progress Reporting Policy (tqdm)

User-visible long-running loops should use `SpliceGrapher.shared.progress.ProgressIndicator`.
This compatibility class now routes through `tqdm` and applies TTY-aware defaults:

- interactive terminal + `verbose=True`: tqdm progress is enabled
- non-interactive stderr or `verbose=False`: progress output is suppressed

Example:

```python
from SpliceGrapher.shared.progress import ProgressIndicator

indicator = ProgressIndicator(1000000, description="loading", verbose=True)
for _ in range(5000):
    indicator.update()
indicator.finish()
```

## Logging Policy (structlog)

Structured logging is canonical in SGN. New/modified runtime modules should use
`structlog` event logging, not stdlib `logging`, `loguru`, or ad-hoc `print`/`stderr`
output for operational messages.

Baseline bootstrap helpers now live in `SpliceGrapher/shared/logging_utils.py`:

- `configure_logging()` for one-time structlog setup
- `get_logger(__name__)` for module-level logger creation

Policy enforcement for touched code is configured in `pyproject.toml` via Ruff
`banned-api` rules.

## Conda / Mamba / Miniconda / Bioconda-Friendly Guidance

Use conda-based environments when needed for scientific stacks:

```bash
mamba create -n splicegrapher-next python=3.12 -y
mamba activate splicegrapher-next
mamba install -c conda-forge -c bioconda numpy scipy matplotlib pysam gffutils networkx -y
python -m pip install -e .
```

Quick compatibility smoke check in that environment:

```bash
python -c "import tqdm; import SpliceGrapher; print(SpliceGrapher.__name__)"
```

`uv` remains the preferred local development workflow in this repository.

## Docker / Nextflow / BioContainers-Friendly Guidance

This repository tracks container/workflow friendliness as a compatibility
requirement during migration. Detailed contract and checklists live in
`docs/workflow-integrations.md`.

Quick container smoke command:

```bash
docker run --rm -v "$PWD":/work -w /work python:3.12-slim bash -lc '
python -m pip install --upgrade pip &&
python -m pip install . &&
python -c "import importlib.metadata as md; print(md.version(\"splicegrapher-next\"))"'
```

## Tutorial Dataset Curation

Canonical modern tutorial/e2e dataset selections (plant + human) are tracked
in:

- `docs/tutorials/dataset_manifest.toml`
- `docs/tutorials/datasets.md`

## Jupyter Tutorial Notebooks

Notebook tutorials and execution notes are in:

- `docs/notebooks/README.md`
- `docs/notebooks/01_legacy_tutorial_foundation.ipynb`
- `docs/notebooks/02_modern_dataset_track.ipynb`

Quick local smoke execution:

```bash
uv sync --group dev
uv run pytest tests/test_notebook_smoke.py -q
```

## Migration Status

This project is in active extraction and stabilization.
Current priority is compatibility, correctness, and tests before deeper
cleanup or re-architecture.
Clean-invariant ratchet policy is enforced by
`scripts/ci/check_clean_invariant.py`:
- Ruff per-file-ignore debt may decrease, but may not grow.
- Executable bits are disallowed on library Python modules unless explicitly
  allowlisted as true scripts/entrypoints.

## License and Provenance

This repository is a continuation/derivative effort and may contain mixed
provenance over time. Preserve original notices and license headers in copied
files. See `LICENSE`, `NOTICE`, `PROVENANCE.md`, and `LICENSES/README.md`.

## Original SpliceGrapher References

Rogers MF, Thomas J, Reddy AS, Ben-Hur A. SpliceGrapher: detecting patterns of
alternative splicing from RNA-Seq data in the context of gene models and EST
data. *Genome Biology*. 2012 Jan 31;13(1):R4. doi:10.1186/gb-2012-13-1-r4.
PMID: 22293517; PMCID: PMC3334585. <https://doi.org/10.1186/gb-2012-13-1-r4>

Rogers MF, Boucher C, Ben-Hur A. SpliceGrapherXT: From Splice Graphs to
Transcripts Using RNA-Seq. In: *Proceedings of the International Conference on
Bioinformatics, Computational Biology and Biomedical Informatics (BCB 2013)*.
Association for Computing Machinery; 2013:247-255.
doi:10.1145/2506583.2506625. <https://doi.org/10.1145/2506583.2506625>
