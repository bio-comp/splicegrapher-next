# splicegrapher-next

`splicegrapher-next` is an unofficial Python 3 modernization/continuation of
SpliceGrapher for shared use by iDiffIR and TAPIS during active migration.

- Repository/distribution name: `splicegrapher-next`
- Python import namespace (compatibility): `SpliceGrapher`

## Why This Repo Exists

Historically, iDiffIR and TAPIS relied on divergent SpliceGrapher copies.
This repository is the compatibility bridge that reduces duplication and
centralizes migration work while preserving behavior and imports.

## Phase 1 Scope

Goals:

- Extract and stabilize modernized SpliceGrapher components from iDiffIR.
- Preserve `SpliceGrapher` import compatibility where practical.
- Keep visualization (`view/`, `plot/`) as first-class functionality.
- Enable incremental external adoption by iDiffIR and TAPIS.

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

### Editable pip Fallback

```bash
python -m pip install -e .
```

### Import Compatibility Check

```bash
uv run python -c "import SpliceGrapher; print(SpliceGrapher.__name__)"
```

## Conda / Mamba / Miniconda / Bioconda-Friendly Guidance

Use conda-based environments when needed for scientific stacks:

```bash
mamba create -n splicegrapher-next python=3.12 -y
mamba activate splicegrapher-next
mamba install -c conda-forge -c bioconda numpy scipy matplotlib pysam gffutils networkx -y
python -m pip install -e .
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

## Migration Status

This project is in active extraction and stabilization.
Current priority is compatibility, correctness, and tests before deeper
cleanup or re-architecture.

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
