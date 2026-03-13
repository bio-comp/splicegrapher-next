# Polars GFF Hard-Clean Design

## Goal
Modernize `SpliceGrapher/formats/polars_gff.py` by removing `object`-typed runtime helper boundaries and tightening the module around explicit typed feature protocols.

## Motivation
- `polars_gff.py` still relies on `object` in runtime helper seams such as feature coordinate extraction and transcript iteration.
- SGN now treats `Any` and `object` in runtime code as code smells unless the value is genuinely unconstrained.
- The module is isolated enough to hard-clean without reopening the larger parser and gene-model cuts that already landed.

## Scope
- remove `object` from runtime helper signatures where a protocol or concrete type can express the contract
- keep optional-Polars behavior and public helpers intact
- preserve current row flattening behavior for gene-model-backed analytics
- add or update tests to pin the tightened helper boundary

## Constraints
- no compatibility wrapper layer
- keep local-only tracker and AGENTS files unstaged
