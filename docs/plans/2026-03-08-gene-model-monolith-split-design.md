# Gene-Model Monolith Split Design

Date: 2026-03-08
Issue: #175
Branch: `sg-next/175-genemodel-monolith-plan`

## Goal

Inventory and plan the split of the two largest remaining SGN gene-model modules:

- `SpliceGrapher/formats/models.py` (1018 lines)
- `SpliceGrapher/formats/parsers/gene_model_gff.py` (903 lines)

This tranche is planning-first. It does not change runtime behavior.

## Current State

### `SpliceGrapher/formats/models.py`

This file mixes several distinct responsibilities:

1. Domain constants and type aliases
2. Interval and sort helpers
3. Per-chromosome interval index logic
4. Genomic coordinate abstraction (`Locus`)
5. Feature entities (`BaseFeature`, `TranscriptRegion`, `Exon`, `CDS`, UTRs)
6. Higher-level containers (`Transcript`, `Gene`, `PseudoGene`)

The module has only a small direct import surface inside SGN:
- `SpliceGrapher/formats/gene_model.py`
- `SpliceGrapher/formats/serializers.py`
- `SpliceGrapher/formats/parsers/gene_model_gff.py`
- `tests/test_gene_model_serializers.py`

That means the real compatibility boundary is already `SpliceGrapher/formats/gene_model.py`, not `models.py` itself.

### `SpliceGrapher/formats/parsers/gene_model_gff.py`

This file is already a parser boundary, but it still bundles multiple concerns:

1. Parser protocols and stats containers
2. annotation/chromosome normalization helpers
3. pending-child queues and parent-cache mechanics
4. parent resolution logic
5. per-record handlers (`gene`, `mrna`, `exon`, CDS/UTR, misc, chromosome)
6. top-level loading orchestration

Its direct runtime import surface is very small:
- `SpliceGrapher/formats/gene_model.py` imports only `load_gene_model_records`
- tests import parser internals directly for bounded-cache behavior

## Design Decision

### 1. Keep `gene_model.py` as the public orchestration façade

Do not spend this tranche breaking the `GeneModel` entry surface.
The existing façade is already the right compatibility boundary for staged cleanup.

### 2. Split `models.py` behind the same module path

Do not hard-delete `SpliceGrapher/formats/models.py` yet.
Instead, convert it into a thin aggregation module over smaller implementation files.

Recommended target structure:

- `SpliceGrapher/formats/model_domain.py`
  - domain constants, aliases, `Locus`
  - feature entities (`BaseFeature`, `TranscriptRegion`, `Exon`, `CDS`, `FpUtr`, `TpUtr`)
  - higher-level entities (`Transcript`, `Gene`, `PseudoGene`)
- `SpliceGrapher/formats/model_index.py`
  - `IntervalQuery`
  - `ChromosomeGeneIndex`
  - `Chromosome`
  - sort/filter/interval helper functions
- `SpliceGrapher/formats/models.py`
  - compatibility aggregation only

Reason:
- lowers blast radius because direct imports of `SpliceGrapher.formats.models` still work
- matches the already-established façade pattern used in `gene_model.py`
- avoids scattering constants into many tiny modules

### 3. Split `gene_model_gff.py` into parser internals plus a thin boundary file

Recommended target structure:

- `SpliceGrapher/formats/parsers/gene_model_gff_context.py`
  - `GeneModelLike`
  - `ParseStats`
  - `ParsedRecord`
  - `ParseContext`
- `SpliceGrapher/formats/parsers/gene_model_gff_records.py`
  - annotation parsing helpers
  - record-line parsing
  - shared record utility functions
- `SpliceGrapher/formats/parsers/gene_model_gff_handlers.py`
  - gene/mrna/exon/transcript-region/misc/chromosome handlers
  - parent resolution helpers and pending-child drainage
- `SpliceGrapher/formats/parsers/gene_model_gff.py`
  - `load_gene_model_records`
  - orchestration only

Reason:
- preserves the current parser import path used by `gene_model.py`
- isolates the handler complexity without forcing a public API break
- keeps bounded-cache internals testable in focused modules

## Compatibility Boundary

- Keep `SpliceGrapher/formats/gene_model.py` stable during the split.
- Keep `SpliceGrapher/formats/parsers/gene_model_gff.py` exporting `load_gene_model_records`.
- Keep `SpliceGrapher/formats/models.py` importable while shrinking it into an aggregator.

## Risks

1. Internal tests currently assert direct source-level dependencies, e.g. parser importing `formats.models` directly.
2. Serializer code imports `formats.models` directly, so partial splits must avoid circular imports.
3. Constants and aliases should not be scattered widely; the split must remain responsibility-based, not file-count-maximizing.

## TDD Strategy

1. Add/import rewrite tests for the future split seams where appropriate.
2. Split `models.py` first while keeping `formats.models` as the stable import surface.
3. Split `gene_model_gff.py` second, keeping `load_gene_model_records` at the existing path.
4. Run targeted gene-model/parser/serializer tests after each step.

## Non-Goals

- No new extraction work
- No `GeneModel` API redesign in this tranche
- No behavior changes in parsing or serialization logic during the first split
