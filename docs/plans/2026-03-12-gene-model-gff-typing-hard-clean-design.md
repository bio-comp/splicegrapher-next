# Gene Model GFF Helper Typing Hard-Clean Design

## Goal
Remove the remaining object-typed helper boundaries from the gene-model GFF parser helper layer.

## Motivation
- `gene_model_gff_records.py` still exposes `annotation_value(..., key: object)`.
- `gene_model_gff_resolution.py` still uses `Mapping[str, object]` helper boundaries.
- These are now among the clearest remaining SGN runtime typing smells after the recent package and parser cleanup work.

## Scope
- replace `object`-typed helper seams with precise types or typed protocols
- preserve parser behavior and existing loading semantics
- add or tighten tests to pin the helper boundary

## Constraints
- no compatibility wrapper layer
- keep local-only tracker and AGENTS files unstaged
