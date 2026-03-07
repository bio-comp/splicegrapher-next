# ADR: SpliceGraph Topology Backend and Boundary

- Date: 2026-03-06
- Status: Accepted
- Issue: `#168`

## Context

`SpliceGrapher/SpliceGraph.py` currently mixes four concerns in one mutable object graph:

1. genomic node data,
2. graph topology storage,
3. traversal/mutation APIs, and
4. GFF serialization.

That design is encoded through `SpliceGraph.nodeDict`, node-local `parents` and `children` lists, and `SpliceGraph.writeGFF()`. The result is high coupling, duplicated traversal logic, a broad mutation surface, and no clean serialization boundary.

A staged compatibility mirror was explicitly rejected for this rewrite. Maintaining a `networkx` graph alongside node-local `parents` / `children` state would create a split-brain architecture with multiple sources of truth.

## Decision

For issue `#168`, `SpliceGraph` topology is owned by `networkx.DiGraph`.

The architecture boundary is:

- `SpliceGraph` is the sole owner of graph topology.
- `SpliceGraphNode` is a pure data container for genomic coordinates, identifiers, strand/chromosome metadata, and node attributes.
- Nodes do not own or mutate parent/child relationships.
- All topology mutations must route through `SpliceGraph` APIs.
- All topology traversals must route through `SpliceGraph` APIs.
- Serialization is not a `SpliceGraph` responsibility; GFF writing moves to `SpliceGrapher/formats/writers/splice_graph.py`.

Concretely, the rewrite should:

- replace topology storage with a `networkx.DiGraph`,
- remove `nodeDict` as the canonical graph store,
- delete node-local `parents` and `children` state,
- delete node-level edge mutation helpers such as `addChild`, `addParent`, `removeChild`, and `removeParent`,
- expose graph-owned APIs for operations such as adding/removing edges, removing nodes, enumerating roots/leaves, and predecessor/successor traversal,
- delete `SpliceGraph.writeGFF()` from the data model, and
- provide a writer boundary in `SpliceGrapher/formats/writers/splice_graph.py` that reuses a shared output helper instead of duplicating writer-local `_open_output` logic.

## Consequences

### Positive

- one source of truth for topology,
- less duplicated traversal code,
- clearer mutation boundaries,
- easier graph algorithms through `networkx`,
- simpler node model semantics, and
- a clean serialization boundary aligned with the existing `formats/writers/` direction.

### Negative

- this is a deliberate breaking change for downstream callers that directly mutate `node.parents` / `node.children`,
- SGN, iDiffIR, and TAPIS call sites will need follow-on migration work,
- `graph.writeGFF(...)` callers will need migration to the new writer API, and
- legacy convenience patterns on nodes will be removed rather than preserved.

## Non-Goals

This ADR does not require parser extraction or a compatibility shim layer.

The rewrite should not preserve node-level topology mutation APIs behind proxies or mirrors. Downstream code should be updated to use graph-owned APIs instead.
