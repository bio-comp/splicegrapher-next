# #185 Design: Core SpliceGraph snake_case hard cut

## Summary

`SpliceGrapher/core/splice_graph.py` still exposes a dense legacy camelCase API even after the networkx rewrite. This slice removes that residue directly by renaming the core graph and node surface to snake_case and updating every SGN caller in the same PR.

## Scope

In scope:
- rename camelCase methods and properties in `SpliceGrapher/core/splice_graph.py`
- update SGN runtime callers, tests, and writers/parsers to the new names
- remove legacy names without compatibility aliases

Out of scope:
- `SpliceGrapher/formats/junction.py` and other non-splice-graph camelCase residue
- broader graph behavior changes
- parser loading architecture beyond the minimal caller updates required for renames

## Constraints

- Keep the current networkx-backed topology model intact.
- Preserve behavior; this is a naming hard cut, not a graph-engine rewrite.
- Treat the parser/writer/test fallout as part of the same slice.
- Do not stage local-only tracker files.

## Rename map

### SpliceGraphNode
- `acceptorEnd` -> `acceptor_end`
- `donorEnd` -> `donor_end`
- `addAltForm` -> `add_alt_form`
- `removeAltForm` -> `remove_alt_form`
- `addAttribute` -> `add_attribute`
- `addCodon` -> `add_codon`
- `addStartCodon` -> `add_start_codon`
- `addEndCodon` -> `add_end_codon`
- `addFormsFromString` -> `add_forms_from_string`
- `addIsoform` -> `add_isoform`
- `addIsoformString` -> `add_isoform_string`
- `altForms` -> `alt_forms`
- `altFormString` -> `alt_form_string`
- `attributeString` -> `attribute_string`
- `codonString` -> `codon_string`
- `downstreamOf` -> `downstream_of`
- `endCodons` -> `end_codons`
- `endCodonString` -> `end_codon_string`
- `hasAS` -> `has_as`
- `hasDisposition` -> `has_disposition`
- `isAltAcceptor` -> `is_alt_acceptor`
- `isAltDonor` -> `is_alt_donor`
- `isKnown` -> `is_known`
- `isPredicted` -> `is_predicted`
- `isRetainedIntron` -> `is_retained_intron`
- `isSkippedExon` -> `is_skipped_exon`
- `isUnresolved` -> `is_unresolved`
- `isoformList` -> `isoform_list`
- `isoformString` -> `isoform_string`
- `putativeChildren` -> `putative_children`
- `putativeParents` -> `putative_parents`
- `startCodons` -> `start_codons`
- `startCodonString` -> `start_codon_string`
- `upstreamOf` -> `upstream_of`
- `altFormSet` -> `alt_form_set`
- `isoformSet` -> `isoform_set`
- `origStart` -> `orig_start`
- `origEnd` -> `orig_end`

### SpliceGraph
- `nodeDict` -> `node_dict`
- `addNode` -> `add_node`
- `addCodons` -> `add_codons`
- `addEdge` -> `add_edge`
- `addEndCodons` -> `add_end_codons`
- `addStartCodons` -> `add_start_codons`
- `attributeString` -> `attribute_string`
- `deleteNode` -> `delete_node`
- `getLeaves` -> `get_leaves`
- `getName` -> `get_name`
- `getNode` -> `get_node`
- `getRoots` -> `get_roots`
- `isEmpty` -> `is_empty`
- `resolvedNodes` -> `resolved_nodes`
- `setName` -> `set_name`
- `unresolvedNodes` -> `unresolved_nodes`

## Caller blast radius

Confirmed callers are concentrated in:
- `SpliceGrapher/formats/parsers/splice_graph.py`
- `SpliceGrapher/formats/writers/splice_graph.py`
- `SpliceGrapher/core/splicing_events.py`
- `SpliceGrapher/core/graph_math.py`
- splice-graph and integration tests

That is a contained hard cut and does not require another package-boundary rethink.
