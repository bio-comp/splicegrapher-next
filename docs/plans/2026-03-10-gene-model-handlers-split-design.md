# Gene-Model Parser Handlers Split Design

Issue: #193  
Branch: `sg-next/193-gene-model-handlers-split`

## Problem

`SpliceGrapher/formats/parsers/gene_model_gff_handlers.py` is the next dense SGN parser monolith. It mixes:
- pending-child queue helpers
- parent / isoform / strand resolution
- record-type dispatch handlers

That keeps the parser facade stable, but the internal file is still too dense to maintain cleanly.

## Goal

Split the handlers module into smaller internal parser modules while preserving the existing `SpliceGrapher.formats.parsers.gene_model_gff` facade.

## Non-Goals

- no parser behavior changes
- no `annotation_io` changes
- no `GeneModel` facade changes

## Target Shape

`SpliceGrapher/formats/parsers/`
- `gene_model_gff_handlers.py` removed
- `gene_model_gff_pending.py`
- `gene_model_gff_resolution.py`
- `gene_model_gff_record_handlers.py`

## Constraints

- keep `load_gene_model_records(...)` behavior stable
- keep parser facade imports stable from `SpliceGrapher.formats.parsers.gene_model_gff`
- use layout tests to pin the new internal structure

## Verification

- focused parser/gene-model test slice
- full `pytest`
- `ruff`, `mypy`, `uv build`
