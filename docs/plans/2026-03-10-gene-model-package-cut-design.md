# GeneModel Package Cut Design

## Goal
Convert `SpliceGrapher/formats/gene_model.py` from a flat facade module into a package boundary while keeping `SpliceGrapher.formats.gene_model` import-stable for SGN callers.

## Motivation
- `gene_model.py` is still the largest remaining SGN runtime facade.
- It mixes three responsibilities:
  - domain symbol re-exports
  - repository/load-write boundary methods
  - `GeneModel` container and query API
- Comparable SGN surfaces have already been cleaned into packages:
  - `formats.models`
  - `formats.fasta`
  - `formats.alignment_io`

## Scope
- Replace `SpliceGrapher/formats/gene_model.py` with `SpliceGrapher/formats/gene_model/`.
- Split internal responsibilities across smaller modules.
- Preserve runtime import behavior for `from SpliceGrapher.formats.gene_model import GeneModel` and existing symbol imports.
- Update layout tests to pin the package boundary.

## Proposed Package Shape
- `SpliceGrapher/formats/gene_model/__init__.py`
  - public facade and re-exports
- `SpliceGrapher/formats/gene_model/repository.py`
  - `GeneModelRepository`
- `SpliceGrapher/formats/gene_model/model.py`
  - `GeneModel`
- optional small helper module only if needed during extraction

## Constraints
- No compatibility shim file left behind at `SpliceGrapher/formats/gene_model.py`.
- No behavior change in gene-model loading/writing/query results.
- Keep local-only tracker and AGENTS files unstaged.
