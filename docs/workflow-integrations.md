# Docker, Nextflow, and BioContainers Friendliness

This document defines the current compatibility contract for containerized and
workflow-orchestrated use of `splicegrapher-next` during extraction/migration.

## Scope

- Keep this repository installable and runnable in non-interactive containers.
- Keep dependency choices friendly to conda-forge/bioconda environments.
- Provide a minimal integration contract for downstream Nextflow wrappers.

This phase is compatibility-first. It does not yet provide a full standalone
pipeline image.

## Docker Friendliness Baseline

Required conventions for new/touched runtime paths:

- Non-interactive execution only (no TTY prompts in library logic).
- Bind-mount friendly paths (`/work` style) for input/output.
- Deterministic dependency setup from `pyproject.toml`.

Reference smoke command (local source mount):

```bash
docker run --rm -v "$PWD":/work -w /work python:3.12-slim bash -lc '
python -m pip install --upgrade pip &&
python -m pip install . &&
python -c "import importlib.metadata as md; print(md.version(\"splicegrapher-next\"))"'
```

## Nextflow Integration Contract (Initial)

Downstream pipelines (iDiffIR/TAPIS) should treat `splicegrapher-next` as a
library dependency with explicit environment setup per process.

Initial process expectations:

- Container or conda environment resolves runtime dependencies.
- Process installs project package from a wheel/sdist or repository checkout.
- Process commands are non-interactive and emit clear stderr diagnostics on
  failure.

Minimal process shape:

```nextflow
process splicegrapher_smoke {
  container 'python:3.12-slim'

  input:
  path repo_root

  script:
  """
  cd ${repo_root}
  python -m pip install .
  python -c "import importlib.metadata as md; print(md.version('splicegrapher-next'))"
  """
}
```

## BioContainers Readiness Checklist

Readiness targets for this repository:

- Runtime dependencies stay available in conda-forge/bioconda ecosystems.
- Avoid unnecessary strict pins that reduce solver compatibility.
- Keep packaging metadata explicit in `pyproject.toml`.
- Provide at least one deterministic smoke command for post-build validation.

Current smoke command target:

- Install package in clean container env.
- Verify distribution metadata resolves (`splicegrapher-next` version lookup).

## Relationship to Other Standards

- Logging standardization is tracked separately in issue `#9`.
- Progress bar/conda standardization is tracked separately in issue `#10`.

This document only covers Docker/Nextflow/BioContainers friendliness for issue
`#11`.
