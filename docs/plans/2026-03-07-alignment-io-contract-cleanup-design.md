# alignment_io.py Cleanup Tranche D Design

## Goal

Harden `SpliceGrapher/formats/alignment_io.py` API clarity and tuple-shape guarantees without breaking the current public ABI.

## Scope

This tranche stays strictly inside `#165`:
- keep public camelCase entrypoints intact
- keep tuple return types intact
- document and test the `alignments=True` contract
- fix depths-file fallback so it preserves the requested tuple shape
- tighten touched typing where it improves local clarity without creating cross-module coupling
- reorganize internal helpers into clearer concern blocks within the same file

Out of scope:
- moving benchmark code
- `SpliceGraph.py` lowercase/shim migration
- introducing snake_case public aliases
- breaking tuple returns into dataclasses or named tuples
- broad extraction or package-layout moves

## Current State

`alignment_io.py` already has a useful internal split:
- source normalization and in-memory SAM wrapping
- pysam-backed alignment opening
- junction/depth/alignment collection
- legacy public wrappers

The primary contract hazard is in `getSamReadData(...)`.

When `alignments=True`, callers reasonably expect a 3-tuple:
- depths map
- junction map
- alignment map

The pysam-backed path satisfies that contract. The depths-file fallback does not; it currently returns only a 2-tuple. That is an ABI footgun because downstream callers can legally unpack three values and fail at runtime only on one branch.

## Chosen Approach

### 1. Preserve public ABI exactly

The module will keep its public camelCase names:
- `getSamReadData`
- `getSamDepths`
- `getSamJunctions`
- `getSamAlignments`
- `getSamHeaders`
- `getSamHeaderInfo`
- `getSamSequences`
- `pysamReadDepths`
- `pysamStrand`
- `isBamFile`
- `isCramFile`
- `makeChromosomeSet`

No new public aliases land in this tranche.

### 2. Make tuple-shape behavior explicit and test-backed

Add coverage for:
- `getSamReadData(..., alignments=True)` on the pysam path
- `getSamReadData(..., alignments=True)` on the depths-file fallback path

The fallback will return `(depths, junctions, {})` when `alignments=True` so both branches satisfy the same structural contract.

### 3. Tighten typing without coupling to the domain model

`pysamReadDepths` should not accept an untyped `gene` parameter. Instead, define a minimal structural protocol describing the attributes it actually uses:
- `id`
- `strand`
- `minpos`
- `maxpos`

This keeps the I/O layer decoupled from `Gene` or `GeneModel` imports.

For source handling, prefer more specific type aliases over raw `object` where practical, but do not force a large rewrite of every helper signature if it makes the file noisier than it makes it safer.

### 4. Reorganize internal code by concern

Within the same file, arrange helpers into four blocks:
1. source/path normalization
2. pysam opening and collection
3. depths-file fallback bridge
4. public legacy API wrappers

This is an internal compartmentalization cut, not a file split.

## Testing Strategy

Regression coverage should prove:
- BAM/SAM path still returns a 3-tuple when `alignments=True`
- alignment map is populated on the pysam path
- depths-file fallback returns a 3-tuple with an empty alignment map when `alignments=True`
- existing 2-tuple behavior remains unchanged when `alignments=False`

Verification remains the usual touched-file gates:
- Ruff check
- Ruff format check
- MyPy on touched modules/tests
- targeted pytest for alignment I/O coverage

## Follow-up Work (Not This Tranche)

1. Move `polars_gff_benchmark.py` out of the runtime package.
2. Introduce a lowercase `splice_graph.py` module with `SpliceGraph.py` as a compatibility shim.
3. Consider future snake_case aliasing or public API cleanup for `alignment_io.py` once downstream callers are ready.
4. Split `alignment_io.py` into smaller modules after the ABI contract is fully pinned by tests.
