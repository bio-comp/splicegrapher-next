"""
Module that encapsulates a splice graph.
"""

from __future__ import annotations

import sys
from collections.abc import Iterable, Iterator
from typing import TextIO

import structlog

from SpliceGrapher.core.enum_coercion import coerce_enum
from SpliceGrapher.core.enums import (
    ALT_SPLICE_EVENT_CODE_BY_NAME,
    ALT_SPLICE_EVENT_NAME_BY_CODE,
    AlternativeSplicingEvent,
    AlternativeSplicingEventName,
    AttrKey,
    EdgeType,
    NodeDisposition,
    RecordType,
    Strand,
)
from SpliceGrapher.core.interval_helpers import (
    InMemoryIntervalIndex,
    interval_contains,
    intervals_overlap,
)
from SpliceGrapher.shared.collection_utils import as_list
from SpliceGrapher.shared.file_utils import ez_open
from SpliceGrapher.shared.format_utils import list_string
from SpliceGrapher.shared.process_utils import idFactory
from SpliceGrapher.shared.progress import ProgressIndicator

LOGGER = structlog.get_logger(__name__)

# Attributes used in standard GFF3 files:
SOURCE_NAME = "SpliceGraph"
DEPRECATED_GENE_REC = RecordType.CLUSTER.value
GENE_REC = RecordType.GRAPH.value
PARENT_REC = RecordType.PARENT.value
CHILD_REC = RecordType.CHILD.value
VALID_GENES = [DEPRECATED_GENE_REC, GENE_REC]
VALID_RECTYPES = [DEPRECATED_GENE_REC, GENE_REC, PARENT_REC, CHILD_REC]

# Graph node attributes for splice graph GFF files
PARENT_ATTR = AttrKey.PARENT.value
ID_ATTR = AttrKey.ID.value
KNOWN_ATTRS = [PARENT_ATTR, ID_ATTR]

# Specific node attribute and values for AS forms:
AS_KEY = AttrKey.ALT_FORM.value

# Attributes for unresolved nodes:
PUTATIVE_PARENTS = "putative_parents"
PUTATIVE_CHILDREN = "putative_children"

# Attributes for predicting protein sequences:
START_CODON_KEY = AttrKey.START_CODON.value
END_CODON_KEY = AttrKey.END_CODON.value

# Alternative splicing event names and abbreviations are defined canonically in core.enums.
IR_NAME = AlternativeSplicingEventName.IR.value
ALT5_NAME = AlternativeSplicingEventName.ALT5.value
ALT3_NAME = AlternativeSplicingEventName.ALT3.value
ALTB3_NAME = AlternativeSplicingEventName.ALTB3.value
ALTB5_NAME = AlternativeSplicingEventName.ALTB5.value
ES_NAME = AlternativeSplicingEventName.ES.value
ALTI_NAME = AlternativeSplicingEventName.ALTI.value
ALTT_NAME = AlternativeSplicingEventName.ALTT.value
AS_NAMES = [event_name.value for event_name in AlternativeSplicingEventName]

IR_ABBREV = AlternativeSplicingEvent.IR.value
ALT5_ABBREV = AlternativeSplicingEvent.ALT5.value
ALT3_ABBREV = AlternativeSplicingEvent.ALT3.value
ALTB3_ABBREV = AlternativeSplicingEvent.ALTB3.value
ALTB5_ABBREV = AlternativeSplicingEvent.ALTB5.value
ES_ABBREV = AlternativeSplicingEvent.ES.value
ALTI_ABBREV = AlternativeSplicingEvent.ALTI.value
ALTT_ABBREV = AlternativeSplicingEvent.ALTT.value
AS_ABBREVS = [event.value for event in AlternativeSplicingEvent]
AS_ABBREV_SET = {event.value for event in AlternativeSplicingEvent}

EVENT_NAME = {
    event.value: ALT_SPLICE_EVENT_NAME_BY_CODE[event].value for event in AlternativeSplicingEvent
}
EVENT_ABBREV = {
    event_name.value: event.value for event_name, event in ALT_SPLICE_EVENT_CODE_BY_NAME.items()
}

ALT_35_NAMES = [ALT3_NAME, ALT5_NAME, ALTB3_NAME, ALTB5_NAME]
ALT_35_ABBREVS = [ALT3_ABBREV, ALT5_ABBREV, ALTB3_ABBREV, ALTB5_ABBREV]
NON_35_NAMES = set(AS_NAMES) - set(ALT_35_NAMES)
NON_35_EVENTS = set(AlternativeSplicingEvent) - {
    AlternativeSplicingEvent.ALT3,
    AlternativeSplicingEvent.ALT5,
    AlternativeSplicingEvent.ALTB3,
    AlternativeSplicingEvent.ALTB5,
}
NON_35_ABBREVS = {event.value for event in NON_35_EVENTS}

# Dispositions for nodes in a graph:
DISPOSITION_KEY = AttrKey.DISPOSITION.value
KNOWN_NODE = NodeDisposition.KNOWN.value
PREDICTED_NODE = NodeDisposition.PREDICTED.value
UNRESOLVED_NODE = NodeDisposition.UNRESOLVED.value

# Edge types for adding edges to new nodes
PARENT_EDGE = EdgeType.PARENT.value
CHILD_EDGE = EdgeType.CHILD.value
ALL_EDGE_TYPES = [PARENT_EDGE, CHILD_EDGE]

# Identifies predicted nodes
PNODE_PREFIX = "pred_"

# Attributes used for displaying splice graphs:
ACCEPTORS_KEY = "acceptors"
DONORS_KEY = "donors"
INFORMATION_KEY = "as_info"

# Keys for retained intron evidence cross-referencing
IR_TP_EVIDENCE = "ir_tp_evidence"
IR_FP_EVIDENCE = "ir_fp_evidence"

# Specific node attribute for isoform names
ISO_KEY = AttrKey.ISOFORMS.value

# Specific graph attribute for expanded gene models
ALT_MODEL_KEY = AttrKey.EXPANDED.value

VALID_STRANDS = [Strand.PLUS.value, Strand.MINUS.value]


# ========================================================================
# Function definitions
def _coerce_alt_splicing_event(
    form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
) -> AlternativeSplicingEvent:
    if isinstance(form, AlternativeSplicingEvent):
        return form
    if isinstance(form, AlternativeSplicingEventName):
        return ALT_SPLICE_EVENT_CODE_BY_NAME[form]

    normalized = form.strip()
    if not normalized:
        raise ValueError("Alternative splicing form cannot be empty")

    try:
        return AlternativeSplicingEvent(normalized)
    except ValueError:
        try:
            return ALT_SPLICE_EVENT_CODE_BY_NAME[AlternativeSplicingEventName(normalized)]
        except ValueError as exc:
            raise ValueError(f"Unknown alternative splicing form: {form}") from exc


def acceptor(node: SpliceGraphNode) -> int:
    """Returns the strand-adjusted position at the start of the acceptor dimer."""
    return node.minpos - 2 if node.strand == "+" else node.maxpos


def containsEdge(node: SpliceGraphNode, edge: Edge) -> bool:
    """Returns true if the Node contains the given Edge; false otherwise.  Note
    that the edge must be completely contained within the node (<,>)."""
    return interval_contains(node, edge, strict=True)


def containsNode(edge: Edge, node: SpliceGraphNode) -> bool:
    """Returns true if the Edge contains the given Node; false otherwise.  Note
    that the node must be completely contained within the edge (<,>)."""
    return interval_contains(edge, node, strict=True)


def childEdges(n: SpliceGraphNode) -> set[Edge]:
    """Returns a set of child edges for a node."""
    return {Edge(n, c) for c in n.children}


def donor(node: SpliceGraphNode) -> int:
    """Returns the strand-adjusted position at the start of the donor dimer."""
    return node.maxpos if node.strand == "+" else node.minpos - 2


def edgeSet(G: SpliceGraph) -> set[Edge]:
    """Returns a complete set of edges found in splice graph G.
    This includes duplicate edges between distinct nodes as found
    in alternate 3'/5' events."""
    result: set[Edge] = set()
    for n in G.nodeDict.values():
        result.update(childEdges(n))
    return result


def intronSet(G: SpliceGraph) -> set[tuple[int, int]]:
    """Returns the set of distinct edges found in a splice graph.
    This does not include any duplicate edges."""
    return {(edge.minpos, edge.maxpos) for edge in edgeSet(G)}


def overlap(a: SpliceGraphNode, b: SpliceGraphNode) -> bool:
    """Returns true if the two nodes overlap and are distinct; false otherwise."""
    return a.id != b.id and intervals_overlap(a, b, inclusive=False)


def overlapsAll(A: SpliceGraphNode, nodeList: Iterable[SpliceGraphNode]) -> bool:
    """Returns true if node A overlaps all nodes in the given node list; false otherwise."""
    for n in nodeList:
        if not overlap(A, n):
            return False
    return True


def parentEdges(n: SpliceGraphNode) -> set[Edge]:
    """Returns a set of parent edges for a node."""
    return {Edge(p, n) for p in n.parents}


# ------------------------------------------------------------------------
# AS annotation functions:
def altSiteEventList(
    nodes: Iterable[SpliceGraphNode],
    siteType: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    verbose: bool = False,
) -> list[set[SpliceGraphNode]]:
    """Returns a list of alt. 3' or alt. 5' events.  Each element consists
    of a set of nodes involved in the same event.  Note: assumes the graph's
    nodes have already been annotated (see detectAltAcceptor and detectAltDonor)."""
    site_event = _coerce_alt_splicing_event(siteType)
    if site_event not in {AlternativeSplicingEvent.ALT3, AlternativeSplicingEvent.ALT5}:
        raise ValueError(f"altNodeList called using non-alternative site type: {siteType}")

    def set_string(s: set[SpliceGraphNode]) -> str:
        return ",".join(n.id for n in s)

    def neighbor_nodes(n: SpliceGraphNode) -> set[SpliceGraphNode]:
        return set(n.parents) if site_event == AlternativeSplicingEvent.ALT3 else set(n.children)

    # Find nodes annotated with the given event
    eventNodes = [n for n in nodes if site_event in n.altFormSet]
    if not eventNodes:
        return []
    if verbose:
        modifier = "parents" if site_event == AlternativeSplicingEvent.ALT3 else "children"
        LOGGER.debug(
            "alt_site_event_list_found_nodes",
            count=len(eventNodes),
            site_type=site_event.value,
        )
        for n in eventNodes:
            LOGGER.debug("alt_site_event_list_node", node=str(n))

    # Sorting nodes ensures that overlapping nodes will
    # not be stored in distinct sets to be merged later
    eventNodes.sort(key=lambda node: (node.minpos, node.maxpos, node.id))
    eventList: list[set[SpliceGraphNode]] = []
    overlap_index = InMemoryIntervalIndex(eventNodes)
    stored: set[SpliceGraphNode] = set()
    for n in eventNodes:
        eset = {n}
        # Grab neighboring nodes (share an edge)
        adj_n = neighbor_nodes(n)
        if verbose:
            LOGGER.debug(
                "alt_site_event_list_neighbors",
                node_id=n.id,
                modifier=modifier,
                neighbors=set_string(adj_n),
            )

        # Find other nodes with the same annotation that overlap the current one
        for m in overlap_index.overlaps(n, inclusive=False):
            if m == n or not overlap(m, n):
                continue
            if verbose:
                LOGGER.debug(
                    "alt_site_event_list_overlap_candidate",
                    node_id=n.id,
                    candidate_id=m.id,
                )
            # Grab other node's neighbors
            adj_m = neighbor_nodes(m)

            # If either node overlaps all of the other's parents (alt3)
            # or children (alt5), they are not part of the same event:
            neighborsContained = overlapsAll(m, adj_n) or overlapsAll(n, adj_m)

            # Otherwise, the nodes are part of the same event
            if not neighborsContained:
                eset.add(m)
            elif verbose:
                LOGGER.debug(
                    "alt_site_event_list_neighbor_containment_excluded",
                    node=str(n),
                    candidate=str(m),
                    modifier=modifier,
                )

        # If the nodes overlapping the current node also overlap
        # distinct other nodes, merge the events into one set
        if verbose:
            LOGGER.debug(
                "alt_site_event_list_event_set",
                node_id=n.id,
                event_set=set_string(eset),
            )
        for e in eventList:
            if e & eset:
                if verbose:
                    LOGGER.debug(
                        "alt_site_event_list_merge_sets",
                        event_set=set_string(eset),
                        existing_set=set_string(e),
                    )
                e.update(eset)
                stored.update(eset)
                if verbose:
                    LOGGER.debug(
                        "alt_site_event_list_merged_set",
                        merged_set=set_string(e),
                    )
                break
            elif verbose:
                LOGGER.debug(
                    "alt_site_event_list_no_intersection",
                    event_set=set_string(eset),
                    existing_set=set_string(e),
                )

        if n not in stored:
            if verbose:
                LOGGER.debug(
                    "alt_site_event_list_store_new_set",
                    event_set=set_string(eset),
                )
            eventList.append(eset)

    if verbose:
        LOGGER.debug(
            "alt_site_event_list_final",
            result=";".join([set_string(e) for e in eventList]),
        )
    return eventList


def detectAltAcceptor(
    n: SpliceGraphNode,
    nodes: list[SpliceGraphNode],
    edges: set[Edge],
) -> None:
    """Looks for evidence that identifies the given node as an alternate acceptor."""
    # Inference rules:
    #  - Establish set of all nodes that overlap the given one and are not graph roots
    #  - Get the set of nodes flanking edges that the node contains (retained introns)
    #  - Get the set of all nodes that overlap parents of the given one (other retained introns)
    #  - Subtract those overlapping parents from the first set
    #  - For those that remain, if any have acceptor sites different from
    #    this one, mark it as an alternate acceptor
    if not n.parents:
        return
    allOverlaps = {o for o in nodes if o.parents and overlap(o, n)}
    containedEdges = {e for e in edges if containsEdge(n, e)}
    flankingNodes = {e.child for e in containedEdges if overlapsAll(n, e.child.parents)}
    parentOverlaps = {o for o in nodes if overlapsAll(o, n.parents)}
    overlapSet = allOverlaps - parentOverlaps - flankingNodes

    for o in overlapSet:
        if o.acceptorEnd() != n.acceptorEnd():
            n.addAltForm(AlternativeSplicingEvent.ALT3)
            return


def detectAltBoth(
    nodes: list[SpliceGraphNode],
    verbose: bool = False,
) -> None:
    """Looks for evidence of simultaneous 3'/5' events (alt. both).  Assumes alt. 3' and
    alt. 5' nodes have already been identified."""
    # Using 5' events as a reference...
    # For each 5' (donor site) event, we look at the nodes to see if any of them share
    # the same child (acceptor site).
    #
    # Examples of alt. both event with 4 nodes:
    #       ||||||||||----------------------|||||||||||||
    #       ||||||||----------------------|||||||||||||
    # Example with 6 nodes:
    #       ||||||||||--------------|||||||||||||
    #       ||||||------------|||||||||||||
    #       |||||||||||||--------------|||||||||||||
    # An alt. both event, and an alt. 5' event:
    #       ||||||||||----------------------|||||||||||||  alt.5' only
    #       ||||||||------------------------|||||||||||
    #       ||||||||||||||----------------------|||||||||||||  alt. both
    #
    # Note that the first two donors have edges to the same acceptor site
    # while the third donor has its own unique acceptor site.  Thus the rule to
    # apply is: any node within an alt 5' group that connects to a unique
    # acceptor that is not shared with any other member of that group becomes
    # an alt. both event.
    #
    # An alt. both event, and an alt. 3' event:
    #       ||||||||||----------------------|||||||||||||  alt.3' only
    #       ||||||||||-------------------|||||||||||
    #       ||||||||||||----------------------|||||||||||||  alt. both
    #
    # Note that the first two acceptors have edges to the same donor site
    # while the third donor has its own unique acceptor site.  Thus the rule to
    # apply is: any node within an alt 3' group that connects to a unique
    # donor that is not shared with any other member of that group becomes
    # an alt. both event.
    alt5Events = altSiteEventList(nodes, AlternativeSplicingEvent.ALT5, verbose=verbose)
    for event in alt5Events:
        for n in event:
            nodeAcceptors = {c.acceptorEnd() for c in n.children}
            otherAcceptors = {c.acceptorEnd() for m in event for c in m.children if m != n}
            shared = nodeAcceptors & otherAcceptors
            if not shared:
                n.removeAltForm(AlternativeSplicingEvent.ALT5)
                n.addAltForm(AlternativeSplicingEvent.ALTB5)

    alt3Events = altSiteEventList(nodes, AlternativeSplicingEvent.ALT3)
    for event in alt3Events:
        for n in event:
            nodeDonors = {c.donorEnd() for c in n.parents}
            otherDonors = {c.donorEnd() for m in event for c in m.parents if m != n}
            shared = nodeDonors & otherDonors
            if not shared:
                n.removeAltForm(AlternativeSplicingEvent.ALT3)
                n.addAltForm(AlternativeSplicingEvent.ALTB3)


def detectAltDonor(
    n: SpliceGraphNode,
    nodes: list[SpliceGraphNode],
    edges: set[Edge],
    verbose: bool = False,
) -> None:
    """Looks for evidence that identifies the given node as an alternate donor."""
    # Inference rules:
    #  - Establish set of all nodes that overlap the given one
    #  - Get the set of nodes flanking edges that the node contains (retained introns)
    #  - Get the set of all nodes that overlap all children of the given one
    #    (other retained introns)
    #  - Subtract those overlapping children from the first set
    #  - For those that remain, if any have acceptor sites different from
    #    this one, mark it as an alternate acceptor
    if not n.children:
        return
    allOverlaps = {o for o in nodes if o.children and overlap(o, n)}
    containedEdges = {e for e in edges if containsEdge(n, e)}
    flankingNodes = {e.parent for e in containedEdges if overlapsAll(n, e.parent.children)}
    childOverlaps = {o for o in nodes if overlapsAll(o, n.children)}
    overlapSet = allOverlaps - childOverlaps - flankingNodes

    for o in overlapSet:
        if o.donorEnd() != n.donorEnd():
            n.addAltForm(AlternativeSplicingEvent.ALT5)
            return


def detectRetainedIntron(n: SpliceGraphNode, edges: set[Edge]) -> None:
    """Looks for evidence in the edge list that identifies the given node as a retained intron."""
    for e in edges:
        if containsEdge(n, e):
            n.addAltForm(AlternativeSplicingEvent.IR)
            return


def detectSkippedExon(n: SpliceGraphNode, edges: set[Edge]) -> None:
    """Looks for evidence in the edge list that identifies the given node as a skipped exon."""
    annotation = AlternativeSplicingEvent.ES
    if not n.parents:
        annotation = AlternativeSplicingEvent.ALTI
    if not n.children:
        annotation = AlternativeSplicingEvent.ALTT
    for e in edges:
        if e.minpos < n.minpos and n.maxpos < e.maxpos:
            n.addAltForm(annotation)
            return


# ------------------------------------------------------------------------
# Graph manipulation functions:
def commonAS(A: SpliceGraph, B: SpliceGraph) -> SpliceGraph:
    """Returns a graph that contains just the nodes and AS
    events found in both A and B.  Uses the fact that
    A^B = AUB - A\\B - B\\A"""
    combined = A.union(B)
    [AminB, BminA] = diffAS(A, B)
    subgraph = graphMinusAS(combined, AminB)
    return graphMinusAS(subgraph, BminA)


def diffAS(A: SpliceGraph, B: SpliceGraph) -> list[SpliceGraph]:
    """Detects alternative splicing differences between two graphs.
    Returns two graphs: one that contains the nodes in A with AS not
    found in B and another that contains nodes in B without events in A."""
    return [graphMinusAS(A, B), graphMinusAS(B, A)]


def equivalentGraphs(A, B):
    """Returns True if SpliceGraphs A and B are equivalent; False otherwise.
    This is distinct from the SpliceGraph.__eq__ method as it ignores
    attributes that may be changed in the prediction process."""
    Anodes = A.nodeDict.values()
    Bnodes = B.nodeDict.values()
    if len(Anodes) != len(Bnodes):
        return False
    for n in Anodes:
        Aedges = childEdges(n)
        try:
            idx = Bnodes.index(n)
            o = Bnodes[idx]
            Bedges = childEdges(o)
            if Aedges != Bedges:
                return False
        except ValueError:
            return False
    return True


def getFirstGraph(f: str | TextIO, *, annotate: bool = False, verbose: bool = False) -> SpliceGraph:
    """Returns just the first splice graph found in a file."""
    try:
        result = next(SpliceGraphParser(f, verbose=verbose))
        if annotate:
            result.annotate()
        return result
    except StopIteration:
        raise ValueError("No graph found in %s" % f)


def graphMinusAS(A: SpliceGraph, B: SpliceGraph) -> SpliceGraph:
    """Detects alternative splicing differences between two graphs.
    Returns a graph that contains just the nodes in A that are
    annotated with AS not found in B.  Note that if A is a subgraph
    of B, this will return None."""
    # Get all nodes with AS annotations.  Note that we do
    # not report nodes missing from one graph or the other.
    nodesA = [n for n in A.resolvedNodes() if n.hasAS()]
    result = SpliceGraph(name=A.getName(), chromosome=A.chromosome, strand=A.strand)
    for a in nodesA:
        aAltSet = set(a.altFormSet)
        if not aAltSet:
            continue

        b = B.getNode(a.start, a.end)
        if b:
            bAltSet = set(b.altFormSet)
            aAltSet = aAltSet - bAltSet
            if not aAltSet:
                continue

        node = result.addNode(a.id, a.minpos, a.maxpos)
        for event in aAltSet:
            node.addAltForm(event)
    return result


def graphSubtract(A: SpliceGraph, B: SpliceGraph, *, resolvedOnly: bool = True) -> SpliceGraph:
    """Returns a graph that represents the nodes and edges in A minus those in B."""
    name = "%s\\%s" % (A.getName(), B.getName())
    result = SpliceGraph(name=name, chromosome=A.chromosome, strand=A.strand)
    # Use full graph position range
    result.minpos = A.minpos
    result.maxpos = A.maxpos

    nodesA = A.resolvedNodes() if resolvedOnly else A.nodeDict.values()
    nodesB = B.resolvedNodes() if resolvedOnly else B.nodeDict.values()
    subsetA = []
    for a in nodesA:
        if a not in nodesB:
            newNode = result.addNode(a.id, a.minpos, a.maxpos)
            newNode.attrs = a.attrs
            subsetA.append(a)

    for a in subsetA:
        for c in a.children:
            if c in subsetA:
                result.addEdge(a.id, c.id)

    return result


def consistencyCoefficients(A, B):
    """Deprecated.  Use recall(A,B) instead."""
    return recall(A, B)


def jaccardCoefficients(A, B):
    """Returns the Jaccard coefficients for the similarity between the
    two graphs' vertex and edge sets, respectively.  (Recall that the
    Jaccard coefficient for a set is |A^B|/|AUB|.)  Note that if both
    A and B are the empty set, the Jaccard coefficient is 1."""
    setA = set(A.resolvedNodes())
    setB = set(B.resolvedNodes())
    AandB = setA.intersection(setB)
    AorB = setA.union(setB)
    nodeCoefficient = float(len(AandB)) / len(AorB) if AorB else 1.0

    edgesA = edgeSet(A)
    edgesB = edgeSet(B)
    AandB = edgesA.intersection(edgesB)
    AorB = edgesA.union(edgesB)
    edgeCoefficient = float(len(AandB)) / len(AorB) if AorB else 1.0

    return nodeCoefficient, edgeCoefficient


def nodeString(nodeList):
    """Returns a string representation of all node ids in a list."""
    return ",".join([n.id for n in nodeList])


def recall(A, B):
    """Returns recall values for the similarity between graphs
    A and B for node and edge sets, respectively.  The recall
    for set A relative to B is |A^B|/|B|.  Note that if A and B
    are both empty, the recall value is 1."""
    nodesB = set(B.resolvedNodes())
    AandB = set(A.resolvedNodes()) & nodesB
    nodeCoefficient = float(len(AandB)) / len(nodesB) if nodesB else 1.0

    edgesB = edgeSet(B)
    AandB = edgeSet(A) & edgesB
    edgeCoefficient = float(len(AandB)) / len(edgesB) if edgesB else 1.0

    return nodeCoefficient, edgeCoefficient


def splitFiles(parser, directory=None):
    """Given a parser on a file that contains multiple models, writes each
    model into its own file prefixed with its graph name.  First removes
    special characters from the name.  For example:
        'AT1G01448|alt 3/5' --> 'AT1G01448_alt_3_5_graph.gff'"""
    import os

    def fixName(s):
        result = s.replace(" ", "_")
        result = result.replace("|", "_")
        result = result.replace("/", "_")
        result = result.replace("\\", "_")
        return result.replace("*", "")

    for g in parser.__iter__():
        outName = "%s_graph.gff" % fixName(g.getName())
        outFile = os.path.join(directory, outName) if directory else outName
        g.writeGFF(outFile)


def updateLeaf(A: SpliceGraph, B: SpliceGraph, *, uniqueLeaf: bool = False) -> None:
    """Method that updates leaf nodes in graph A with nodes from graph B prior to
    merging the two graphs.  We look for longer nodes with the same acceptor site."""
    mergeDict = {}
    for leaf in A.getLeaves():
        longer = [
            n
            for n in B.nodeDict.values()
            if n.acceptorEnd() == leaf.acceptorEnd() and len(n) > len(leaf)
        ]
        if not longer:
            continue
        if uniqueLeaf and len(longer) > 1:
            continue
        maxLen = max([len(n) for n in longer])
        longest = [n for n in longer if len(n) == maxLen][0]

        # If the new length matches any other nodes,
        # merge them after the loop has finished
        other = A.getNode(longest.minpos, longest.maxpos)
        if other:
            mergeDict[leaf] = other
        else:
            leaf.update(longest.minpos, longest.maxpos)

    # Merge and remove nodes that would have been duplicates
    for node in mergeDict:
        other = mergeDict[node]
        for p in node.parents:
            A.addEdge(p.id, other.id)
            p.removeChild(node)
        del A.nodeDict[node.id]


def updateRoot(A: SpliceGraph, B: SpliceGraph, *, uniqueRoot: bool = False) -> None:
    """Method that updates root nodes in graph A using nodes from graph B prior to
    merging the two graphs.  We look for longer nodes that have the same donor site."""
    mergeDict = {}
    for root in A.getRoots():
        longer = [
            n for n in B.nodeDict.values() if n.donorEnd() == root.donorEnd() and len(n) > len(root)
        ]
        if not longer:
            continue
        if uniqueRoot and len(longer) > 1:
            continue
        maxLen = max([len(n) for n in longer])
        longest = [n for n in longer if len(n) == maxLen][0]

        # If the new length matches any other nodes,
        # merge them after the loop has finished
        other = A.getNode(longest.minpos, longest.maxpos)
        if other:
            mergeDict[root] = other
        else:
            root.update(longest.minpos, longest.maxpos)

    # Merge and remove nodes that would have been duplicates
    for node in mergeDict:
        other = mergeDict[node]
        for c in node.children:
            A.addEdge(other.id, c.id)
            c.removeParent(node)
        del A.nodeDict[node.id]


# ========================================================================
# Class definitions
class Edge:
    """Encapsulates an edge in the graph.  These are not stored with a splice
    graph, but used for splice graph creation and annotation."""

    def __init__(self, parent: SpliceGraphNode, child: SpliceGraphNode) -> None:
        self.parent = parent
        self.child = child
        self.pos = sorted([parent.minpos, parent.maxpos, child.minpos, child.maxpos])
        # minpos/maxpos are related to the edge itself, not the outer positions
        self.minpos = self.pos[1]
        self.maxpos = self.pos[2]

    def __eq__(self, o: object) -> bool:
        # Using minpos/maxpos allows an edge to be compared with an exon.
        if not isinstance(o, Edge):
            return NotImplemented
        return (
            self.minpos == o.minpos
            and self.maxpos == o.maxpos
            and self.pos[0] == o.pos[0]
            and self.pos[3] == o.pos[3]
        )

    def __lt__(self, o: Edge) -> bool:
        return tuple(self.pos) < tuple(o.pos)

    def __getitem__(self, i: int) -> int:
        return self.pos[i]

    def __hash__(self) -> int:
        return hash(str(self))

    def __len__(self) -> int:
        return self.maxpos - self.minpos

    def overlaps(self, o: Edge) -> bool:
        """Returns true if the edge portions overlap; false otherwise."""
        return self.maxpos > o.minpos and self.minpos < o.maxpos

    def sameEdge(self, o: Edge) -> bool:
        """Returns true if the edge portions are the same; false otherwise."""
        return self.minpos == o.minpos and self.maxpos == o.maxpos

    def __repr__(self) -> str:
        return f"{self.pos[0]},{self.pos[1]},{self.pos[2]},{self.pos[3]}"

    def __str__(self) -> str:
        return f"{self.pos[1]},{self.pos[2]}"


class NullNode:
    """Null node object encapsulates the most basic information about a node."""

    def __init__(self, start: int, end: int) -> None:
        self.start = start
        self.end = end
        self.minpos = min(start, end)
        self.maxpos = max(start, end)


class SpliceGraphNode:
    """This is the node class to use for constructing splice graphs for GFF input/output."""

    def __init__(
        self,
        id: str,
        start: int,
        end: int,
        strand: str | Strand,
        chrom: str,
        parents: list[SpliceGraphNode] | None = None,
        children: list[SpliceGraphNode] | None = None,
    ) -> None:
        self.id = id
        self.strand = coerce_enum(strand, Strand, field="strand").value
        self.chromosome = chrom
        self.minpos = min(start, end)
        self.maxpos = max(start, end)
        (self.start, self.end) = (
            (self.minpos, self.maxpos) if strand == "+" else (self.maxpos, self.minpos)
        )
        self.parents: list[SpliceGraphNode] = list(parents) if parents is not None else []
        self.children: list[SpliceGraphNode] = list(children) if children is not None else []
        self.attrs: dict[str, str | set[int]] = {}
        self.altFormSet: set[AlternativeSplicingEvent] = set()
        self.isoformSet: set[str] = set()
        # Added for adjustable ranges
        self.origStart = self.start
        self.origEnd = self.end

    def acceptorEnd(self):
        """Returns the position of the exon's 5' (upstream) end."""
        return self.start

    def _sync_alt_form_attr(self) -> None:
        self.attrs[AS_KEY] = ",".join(form.value for form in self.altFormSet)

    def addAltForm(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        """Adds an AS form to the node's list of forms."""
        if isinstance(form, str) and not form.strip():
            return
        event = _coerce_alt_splicing_event(form)
        self.altFormSet.add(event)
        for x in self.altFormSet:
            assert len(x.value) > 0
        self._sync_alt_form_attr()

    def addAttribute(self, key, value):
        """Adds the given key-value pair to a node's attribute list.  If the
        key is already in the attribute dictionary, the new value overwrites the old."""
        self.attrs[key] = value

    def addChild(self, c):
        """Adds a child node (edge) to a node.  If the child is already known the child
        list is unchanged."""
        if c not in self.children:
            self.children.append(c)
            c.addParent(self)

    def addCodon(self, codon, codonType):
        """If the codon falls within the node, its start position is added to the set;
        otherwise the list is unchanged."""

        if (not isinstance(codon, tuple)) or (len(codon) != 2):
            raise ValueError(
                "Codons must be provided as a duple of (start,end) positions: received %s"
                % str(codon)
            )

        pos = min(codon) if self.strand == "+" else max(codon)
        if self.contains(pos):
            assert isinstance(pos, int)
            self.attrs.setdefault(codonType, set())
            self.attrs[codonType].add(pos)

    def addEndCodon(self, codon):
        """Adds an end codon to the node."""
        self.addCodon(codon, END_CODON_KEY)

    def addFormsFromString(self, s):
        """Adds AS forms to node list from a string representation."""
        for form in s.split(","):
            self.addAltForm(form.strip())

    def addIsoformString(self, isoformString):
        """Adds a set of isoform names to the node's set of isoforms."""
        for iso in isoformString.split(","):
            self.addIsoform(iso)

    def addIsoform(self, isoform):
        """Adds an isoform name to the node's set of isoforms."""
        if isoform is None:
            raise ValueError("Received illegal isoform for %s" % self.id)
        self.isoformSet.add(isoform)
        self.attrs[ISO_KEY] = ",".join(self.isoformSet)

    def addParent(self, p):
        """Adds a parent node (edge) to a node.  If the parent is already known the parent
        list is unchanged."""
        if p not in self.parents:
            self.parents.append(p)

    def addStartCodon(self, codon):
        """If the codon position falls within the node it is added to the set;
        otherwise the list is unchanged."""
        self.addCodon(codon, START_CODON_KEY)

    def altForms(self):
        """Returns a list of AS forms associated with the node."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return [form.value for form in self.altFormSet]

    def altFormString(self):
        """Returns a string of AS forms associated with the node."""
        try:
            return self.attrs[AS_KEY]
        except KeyError:
            return ""

    def attributeString(self):
        """Returns a string of all attributes associated with this node, suitable for a
        GFF attributes field.  For example:
                  ID=ei_27;Parent=ei_26;AltForm=Alt. 3',Retained Intron"""
        result = ""
        for k in self.attrs.keys():
            if k == AS_KEY and not self.altFormSet:
                continue
            if k == ISO_KEY and not self.isoformSet:
                continue

            if result:
                result += ";"
            if k in [START_CODON_KEY, END_CODON_KEY]:
                result += "%s=%s" % (k, self.codonString(k))
            else:
                result += "%s=%s" % (k, self.attrs[k])
        return result

    def branchingFactor(self):
        return max(len(self.parents), len(self.children))

    def __lt__(self, other: SpliceGraphNode) -> bool:
        """Permits sorting based on minimum node position; ties use max position."""
        if self.minpos == other.minpos:
            return self.maxpos < other.maxpos
        return self.minpos < other.minpos

    def codons(self, codonType):
        """Returns a list of codon positions within the node, or
        an empty list if there are no codons of the given type."""
        try:
            return sorted(self.attrs[codonType])
        except KeyError:
            return []

    def codonString(self, codonType):
        """Returns a list of codon positions within the node as a string,
        or None if there are no codons of the given type."""
        try:
            codons = sorted(self.attrs[codonType])
            return ",".join(["%d" % x for x in codons])
        except KeyError:
            pass

    def contains(self, pos):
        """Returns true if the given position falls within the node; false otherwise."""
        return self.minpos <= pos <= self.maxpos

    def donorEnd(self):
        """Returns the position of the exon's 3' (downstream) end."""
        return self.end

    def downstreamOf(self, pos):
        """Returns true if the node is downstream of the position; false otherwise."""
        return (self.minpos > pos) if self.strand == "+" else (self.maxpos < pos)

    def __eq__(self, other: object) -> bool:
        """Node equality is based solely on its start/end positions."""
        if not isinstance(other, (SpliceGraphNode, NullNode)):
            return NotImplemented
        return self.minpos == other.minpos and self.maxpos == other.maxpos

    def endCodons(self):
        """Returns a list of end codon start positions within the node."""
        return self.codons(END_CODON_KEY)

    def endCodonString(self):
        """Returns a list of end codon start positions within the node as a string."""
        return self.codonString(END_CODON_KEY)

    def gffString(self):
        """Returns a GFF-formatted string representing the given graph feature."""
        if self.parents:
            # Example: chr1	SpliceGraph	child	37373	37398	.	-	.	ID=gm_20;Parent=gm_19
            result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s;Parent=%s" % (
                self.chromosome,
                SOURCE_NAME,
                CHILD_REC,
                self.minpos,
                self.maxpos,
                self.strand,
                self.id,
                nodeString(self.parents),
            )
        else:
            # Example: chr1	SpliceGraph	parent	37569	37757	.	-	.	ID=gm_19
            result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\tID=%s" % (
                self.chromosome,
                SOURCE_NAME,
                PARENT_REC,
                self.minpos,
                self.maxpos,
                self.strand,
                self.id,
            )
        if self.attrs:
            result += ";%s" % self.attributeString()
        return result

    def hasAS(self):
        """Returns true if the node has evidence of AS; false otherwise."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return bool(self.altFormSet)

    def hasDisposition(self, disposition):
        """Returns true if the node has the given disposition; false otherwise."""
        try:
            return self.attrs[DISPOSITION_KEY] == disposition
        except KeyError:
            return False

    def __hash__(self):
        hashString = "%d,%d" % (self.minpos, self.maxpos)
        return hashString.__hash__()

    def isAltAcceptor(self):
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return AlternativeSplicingEvent.ALT3 in self.altFormSet

    def isAltDonor(self):
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return AlternativeSplicingEvent.ALT5 in self.altFormSet

    def isKnown(self):
        """Returns true if the node is known from the gene model; false otherwise."""
        return self.hasDisposition(KNOWN_NODE)

    def isLeaf(self):
        """Returns true if the node is a leaf node; false otherwise."""
        return len(self.children) == 0

    def isPredicted(self):
        """Returns true if the node is predicted; false otherwise."""
        return self.hasDisposition(PREDICTED_NODE)

    def isRetainedIntron(self):
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return AlternativeSplicingEvent.IR in self.altFormSet

    def isRoot(self):
        """Returns true if the node is a root node; false otherwise."""
        return len(self.parents) == 0

    def isSkippedExon(self):
        """Returns true if the node represents a retained intron; false otherwise."""
        for x in self.altFormSet:
            assert len(x.value) > 0
        return AlternativeSplicingEvent.ES in self.altFormSet

    def isUnresolved(self):
        """Returns true if the node is unresolved; false otherwise."""
        return self.hasDisposition(UNRESOLVED_NODE)

    def isoformList(self):
        """Returns a list of isoforms associated with the node."""
        return list(self.isoformSet)

    def isoformString(self):
        """Returns a string representation of isoforms associated with the node."""
        try:
            return self.attrs[ISO_KEY]
        except KeyError:
            return None

    def __len__(self):
        """Returns the length of the genomic region represented by the node."""
        return self.maxpos - self.minpos + 1

    def putativeChildren(self):
        """Returns a set of putative children for an unresolved node,
        or an empty set if there are none or the node is not unresolved."""
        try:
            return self.attrs[PUTATIVE_CHILDREN].split(",")
        except KeyError:
            return set()

    def putativeParents(self):
        """Returns a set of putative parents for an unresolved node,
        or None if there are none or the node is not unresolved."""
        try:
            return self.attrs[PUTATIVE_PARENTS].split(",")
        except KeyError:
            return set()

    def removeAltForm(
        self,
        form: str | AlternativeSplicingEvent | AlternativeSplicingEventName,
    ) -> None:
        """Adds an AS form to the node's list of forms."""
        if isinstance(form, str) and not form.strip():
            return
        event = _coerce_alt_splicing_event(form)
        self.altFormSet.discard(event)
        for x in self.altFormSet:
            assert len(x.value) > 0
        self._sync_alt_form_attr()

    def removeChild(self, child):
        """Removes the child from this node's list."""
        self.children.remove(child)

    def removeParent(self, parent):
        """Removes the parent from this node's list."""
        self.parents.remove(parent)

    def __repr__(self):
        if self.parents:
            return "%s %d-%d children=%s parents=%s" % (
                self.id,
                self.start,
                self.end,
                nodeString(self.children),
                nodeString(self.parents),
            )
        else:
            return "Root %s %d-%d children=%s" % (
                self.id,
                self.start,
                self.end,
                nodeString(self.children),
            )

    def setUnresolved(
        self,
        preserve: list[str] | None = None,
        acceptors: list[int] | None = None,
        donors: list[int] | None = None,
    ) -> None:
        """
        Makes a node unresolved by setting its disposition attribute and
        by removing edges between the node and parents and children in the graph.
        Caller may preserve some display edges by listing nodes in a preserve list.
        """
        preserve_nodes = preserve or []
        acceptor_sites = acceptors or []
        donor_sites = donors or []

        self.addAttribute(DISPOSITION_KEY, UNRESOLVED_NODE)
        if acceptor_sites:
            self.addAttribute(ACCEPTORS_KEY, list_string(acceptor_sites))
        if donor_sites:
            self.addAttribute(DONORS_KEY, list_string(donor_sites))

        for c in list(self.children):
            if c.id in preserve_nodes:
                continue
            c.removeParent(self)
            self.children.remove(c)

        for p in list(self.parents):
            if p.id in preserve_nodes:
                continue
            p.removeChild(self)
            self.parents.remove(p)

    def startCodons(self):
        """Returns a list of start codon start positions within the node."""
        return self.codons(START_CODON_KEY)

    def startCodonString(self):
        """Returns a list of start codon start positions within the node as a string."""
        return self.codonString(START_CODON_KEY)

    def __str__(self):
        if self.parents:
            return "%s %d-%d (parents:%d/children:%d)" % (
                self.id,
                self.start,
                self.end,
                len(self.parents),
                len(self.children),
            )
        else:
            return "Root %s %d-%d (children:%d)" % (
                self.id,
                self.start,
                self.end,
                len(self.children),
            )

    def update(self, minpos, maxpos):
        """Mutates the node based on the revised min/max positions."""
        assert minpos <= maxpos
        self.minpos = minpos
        self.maxpos = maxpos
        self.start = minpos if self.strand == "+" else maxpos
        self.end = maxpos if self.strand == "+" else minpos

    def upstreamOf(self, pos):
        """Returns true if the node is upstream of the position; false otherwise."""
        return (self.maxpos < pos) if self.strand == "+" else (self.minpos > pos)


class SpliceGraph:
    """Main class for storing splice graph data."""

    def __init__(self, name, chromosome, strand):
        """
        Instantiate a splice graph.

        :Parameters:
           'name'       - gene name to associate with the graph
           'chromosome' - chromosome/scaffold name to associate with the
                          graph (GFF format requirement)
           'strand'     - strand associated with the graph
        """
        self.chromosome = chromosome
        self.strand = coerce_enum(strand, Strand, field="strand").value
        self.nodeDict = {}
        self.attrs = {}
        self.minpos = sys.maxsize
        self.maxpos = 0
        # NB: guarantees attributes exist for the graph
        self.setName(name)

    def addNode(self, newId, start, end):
        """Add a node to the splice graph if it doesn't already exist."""
        try:
            # Look for another node with the same start/end positions
            tmpNode = NullNode(start, end)
            allNodes = list(self.nodeDict.values())
            idx = allNodes.index(tmpNode)
            node = allNodes[idx]
        except ValueError:
            self.nodeDict[newId] = SpliceGraphNode(newId, start, end, self.strand, self.chromosome)
            node = self.nodeDict[newId]
            self.minpos = min(self.minpos, node.minpos)
            self.maxpos = max(self.maxpos, node.maxpos)

        return node

    def addCodons(self, codonList, codonType):
        """Adds codons to every node in the graph."""
        for codon in codonList:
            for node in self.nodeDict.values():
                node.addCodon(codon, codonType)

    def addEdge(self, pid, cid):
        """Add an edge to the splice graph if it doesn't already exist."""
        try:
            parent = self.nodeDict[pid]
        except KeyError:
            raise Exception("Error adding edge from node %s: node not found in graph" % pid)

        try:
            child = self.nodeDict[cid]
        except KeyError:
            raise Exception("Error adding edge to node %s: node not found in graph" % cid)

        parent.addChild(child)

    def addEndCodons(self, codonList):
        """Adds end codons to all nodes in the graph."""
        self.addCodons(codonList, END_CODON_KEY)

    def addStartCodons(self, codonList):
        """Adds start codons to all nodes in the graph."""
        self.addCodons(codonList, START_CODON_KEY)

    def adjust(self, adjustment):
        """
        Shifts the positions of all nodes in the graph by the given amount.
        This is useful for example for converting a graph based on positions
        [0,n-1] for an environment that uses positions [1,n].
        """
        for n in self.nodeDict.values():
            n.start += adjustment
            n.end += adjustment
            n.minpos += adjustment
            n.maxpos += adjustment

        self.minpos += adjustment
        self.maxpos += adjustment

    def altForms(self):
        """Returns a set of all AS forms found in the graph."""
        result: set[str] = set()
        for n in self.nodeDict.values():
            result.update(n.altForms())
        for x in result:
            assert len(x) > 0
        return result

    def annotate(self, verbose=False):
        """Updates all AS annotations for the graph."""
        edges = edgeSet(self)
        nodes = self.resolvedNodes()
        for n in nodes:
            # First reset the node AS annotations
            n.altFormSet = set()
            n.attrs[AS_KEY] = ""
            for x in n.altFormSet:
                assert len(x.value) > 0

            # Assign annotations to the node
            detectSkippedExon(n, edges)
            detectRetainedIntron(n, edges)
            if n.children:
                detectAltDonor(n, nodes, edges, verbose=verbose)
            if n.parents:
                detectAltAcceptor(n, nodes, edges)

        # Removed (11/21/2011) after conversation with Asa:
        ## detectAltBoth(nodes)

    def attributeString(self):
        """Returns a string of all graph attributes, for GFF attributes fields."""
        return ";".join(["%s=%s" % (k, self.attrs[k]) for k in sorted(self.attrs.keys())])

    def branchingStats(self):
        """Returns the min, max and average branching factor for the nodes in the graph."""
        branches = [n.branchingFactor() for n in self.resolvedNodes()]
        # NB: len(branches) has to be at least 1 since there must be at least 1 node to have a graph
        try:
            avg = float(sum(branches)) / len(branches)
        except ZeroDivisionError:
            raise ValueError("Attempted to compute branching factor on an empty graph")

        return min(branches), max(branches), avg

    def __lt__(self, other: SpliceGraph) -> bool:
        """Permits sorting graphs based on minimum position; ties use max position."""
        if self.minpos == other.minpos:
            return self.maxpos < other.maxpos
        return self.minpos < other.minpos

    def deleteNode(self, n):
        """Removes a node from a graph, along with all edges attached to it.
        Returns the node that was deleted."""
        # Caller must catch KeyError if this fails
        node = self.nodeDict[n.id]
        for p in node.parents:
            p.removeChild(node)
        for c in node.children:
            c.removeParent(node)
        del self.nodeDict[n.id]
        return node

    def distinctAltEvents(self):
        """Returns a list of distinct AS events found in the graph.
        Each event is reported as a tuple: (node-start,node-end,event-type,count).
        Alt. 5'/Alt. 3' events are the only ones with a count greater than 1 and
        are reported using the shortest node involved."""
        result = []
        alt5Nodes = [
            n for n in self.resolvedNodes() if AlternativeSplicingEvent.ALT5 in n.altFormSet
        ]
        while alt5Nodes:
            n = alt5Nodes.pop()
            others = [o for o in alt5Nodes if o.acceptorEnd() == n.acceptorEnd()]
            for o in others:
                if len(o) < len(n):
                    n = o
                alt5Nodes.remove(o)
            result.append(
                (
                    n.start,
                    n.end,
                    AlternativeSplicingEvent.ALT5.value,
                    len(others) + 1,
                )
            )

        alt3Nodes = [
            n for n in self.resolvedNodes() if AlternativeSplicingEvent.ALT3 in n.altFormSet
        ]
        while alt3Nodes:
            n = alt3Nodes.pop()
            others = [o for o in alt3Nodes if o.donorEnd() == n.donorEnd()]
            for o in others:
                if len(o) < len(n):
                    n = o
                alt3Nodes.remove(o)
            result.append(
                (
                    n.start,
                    n.end,
                    AlternativeSplicingEvent.ALT3.value,
                    len(others) + 1,
                )
            )

        for n in self.resolvedNodes():
            for as_event in NON_35_EVENTS:
                if as_event in n.altFormSet:
                    result.append((n.start, n.end, as_event.value, 1))

        return result

    def downstreamOf(self, a, b):
        """Returns true if a is downstream of b; false otherwise."""
        return (a > b) if self.strand == "+" else (b > a)

    def duplicate(self):
        """Returns a duplicate of the current graph.  This performs
        the same function as copy.deepcopy(), but without the risk
        of stack overflow for large graphs."""
        result = SpliceGraph(self.getName(), self.chromosome, self.strand)
        result.minpos = self.minpos
        result.maxpos = self.maxpos
        result.attrs.update(self.attrs)

        # First pass: create nodes in graph
        for n in self.nodeDict.values():
            newNode = result.addNode(n.id, n.start, n.end)
            newNode.attrs.update(n.attrs)
            newNode.isoformSet = set(n.isoformSet)
            newNode.attrs[ISO_KEY] = ",".join(newNode.isoformSet)

        # Second pass: create edges in graph
        for n in self.nodeDict.values():
            for c in n.children:
                result.addEdge(n.id, c.id)

        return result

    def __eq__(self, other):
        """Returns true if two splice graphs are identical; false otherwise."""
        nodes = self.nodeDict.values()
        otherNodes = other.nodeDict.values()
        if len(nodes) != len(otherNodes):
            return False
        for n in nodes:
            edges = childEdges(n)
            try:
                idx = otherNodes.index(n)
                o = otherNodes[idx]
                otherEdges = childEdges(o)
                if edges != otherEdges:
                    return False
                if n.attrs != o.attrs:
                    return False
            except ValueError:
                return False
        return True

    def expandedModel(self, name, oldRange, newRange):
        """Marks a graph when it expands a gene model beyond its original bounds.
        oldRange and newRange must be (minpos,maxpos) pairs."""
        if len(oldRange) != 2 or not all(isinstance(v, int) for v in oldRange):
            raise ValueError("old range must be 2 int values")

        if len(newRange) != 2 or not all(isinstance(v, int) for v in newRange):
            raise ValueError("new range must be 2 int values")

        self.attrs[ALT_MODEL_KEY] = "%s from (%d,%d) to (%d,%d)" % (
            name,
            oldRange[0],
            oldRange[1],
            newRange[0],
            newRange[1],
        )

    def getAcceptors(self):
        """Returns a list of distinct acceptor sites in the graph."""
        result = {acceptor(n) for n in self.nodeDict.values() if not n.isRoot()}
        return list(result)

    def getDonors(self):
        """Returns a list of distinct donor sites in the graph."""
        result = {donor(n) for n in self.nodeDict.values() if not n.isLeaf()}
        return list(result)

    def getNode(self, start, end):
        """Returns the node represented by the given positions, if it exists."""
        # NB: May want to add option to locate only resolved nodes
        try:
            # Look for a node with the same start/end positions
            tmpNode = NullNode(start, end)
            allNodes = list(self.nodeDict.values())
            idx = allNodes.index(tmpNode)
            return allNodes[idx]
        except ValueError:
            return None

    def getLeaves(self):
        """Returns a list of nodes that have no children."""
        return [n for n in self.nodeDict.values() if not n.children]

    def getName(self):
        """Central method for retrieving the name.  This is a temporary
        hack until a better method is implemented."""
        return self.attrs[ID_ATTR]

    def getRoots(self):
        """Returns a list of nodes that have no parents."""
        return [n for n in self.nodeDict.values() if not n.parents]

    def gffString(self, node=None):
        """Returns a GFF-formatted string representing the given graph feature."""
        result = "%s\t%s\t%s\t%d\t%d\t.\t%s\t.\t%s" % (
            self.chromosome,
            SOURCE_NAME,
            GENE_REC,
            self.minpos,
            self.maxpos,
            self.strand,
            self.attributeString(),
        )
        return result

    def hasAS(self):
        """Returns true if the graph has evidence of AS; false otherwise."""
        for n in self.resolvedNodes():
            if n.hasAS():
                return True
        return False

    def isEmpty(self):
        """Returns true if the graph is empty (has no nodes or edges); false otherwise."""
        return len(self.nodeDict) == 0

    def isoformDict(self):
        """Returns a dictionary where each isoform name is
        keyed to a sorted list of nodes in the isoform."""
        result = {}
        nodes = self.resolvedNodes()
        nodes.sort(
            key=lambda node: (node.minpos, node.maxpos, node.id),
            reverse=(self.strand == "-"),
        )
        for n in nodes:
            for iso in n.isoformSet:
                result.setdefault(iso, [])
                result[iso].append(n)
        return result

    def isoformGraph(self, isoNames, *, fuzzyMatch: bool = False):
        """Returns a splice graph that represents just the given isoforms.
        Isoforms may be given as a comma-separated string, a list or a set.
        Setting 'fuzzyMatch' to True allows you to specify part of an isoform
        name in the list; it will match all isoforms that contain that part."""
        isoNames = as_list(isoNames)
        result = None
        for name in isoNames:
            if fuzzyMatch:
                nodes = []
                query = name.upper()
                for n in self.resolvedNodes():
                    for x in n.isoformSet:
                        if x.upper().find(query) >= 0:
                            nodes.append(n)
            else:
                nodes = [n for n in self.resolvedNodes() if name in n.isoformSet]

            if not nodes:
                raise ValueError(
                    "%s not found in graph %s (isoforms %s)\n"
                    % (name, self.getName(), ",".join(self.isoforms()))
                )

            graph = SpliceGraph(name, self.chromosome, self.strand)
            prev = None
            for n in nodes:
                curr = graph.addNode(n.id, n.minpos, n.maxpos)
                curr.attrs = dict(n.attrs)
                curr.isoformSet = {name}
                curr.attrs[ISO_KEY] = ",".join(curr.isoformSet)
                if prev:
                    graph.addEdge(prev.id, curr.id)
                prev = curr
            result = graph if result is None else result.union(graph)
        if result:
            result.setName("|".join(isoNames))
        return result

    def isoforms(self):
        """Returns a list of all isoform strings associated with nodes in the graph."""
        isoSet = set()
        for n in self.resolvedNodes():
            isoSet = isoSet | n.isoformSet
        return sorted(list(isoSet))

    def __len__(self):
        return self.maxpos - self.minpos + 1 if self.maxpos > self.minpos else 0

    def mergeNodes(self, a, b):
        """Merges two nodes by using the lowest min position and the highest
        max position of the two nodes.  The new node will have the same identifier
        as the first original node.  All parents and all children are combined and the
        original nodes will be deleted from the graph.  Returns the original nodes."""
        node1 = self.nodeDict[a]
        node2 = self.nodeDict[b]
        minpos = min(node1.minpos, node2.minpos)
        maxpos = max(node1.maxpos, node2.maxpos)

        # Caller must catch KeyError if this fails
        self.deleteNode(node1)
        self.deleteNode(node2)

        # Either create a new node or grab an existing one
        newNode = self.addNode(node1.id, minpos, maxpos)

        # Only add parents/children that are still in the graph
        for c in node1.children + node2.children:
            if c.id in self.nodeDict:
                newNode.addChild(c)
        for p in node1.parents + node2.parents:
            if p.id in self.nodeDict:
                p.addChild(newNode)

        self.nodeDict[newNode.id] = newNode
        return (node1, node2)

    def nameString(self):
        """Returns the graph name.  Graphs that have been merged will have more than
        one name given as a hyphenated list."""
        names = set(self.getName().split(","))
        if len(names) > 1:
            return "Merged " + "+".join(names)
        else:
            return self.getName()

    def __ne__(self, other):
        """Returns true if two splice graphs are different; false otherwise."""
        return not self.__eq__(other)

    def __str__(self):
        return "%s (%s) %d-%d (%d nodes)" % (
            self.getName(),
            self.strand,
            self.minpos,
            self.maxpos,
            len(self.nodeDict),
        )

    def resolvedNodes(self):
        """Returns a list of all resolved nodes: those actually used in the graph."""
        return [x for x in self.nodeDict.values() if not x.isUnresolved()]

    def setName(self, name):
        """Provides a central method for storing the name in both places.
        This is a temporary hack until a better method is implemented."""
        self.name = name
        self.attrs[ID_ATTR] = name

    def union(
        self,
        other: SpliceGraph,
        *,
        keepName: bool = False,
        mergeEnds: bool = False,
        idGenerator: Iterator[str] | None = None,
    ) -> SpliceGraph:
        """Merges two graphs into one.  Adds nodes that are distinct and merges nodes
        that have the same start/end positions.  Merged nodes will share AS
        information.  New attributes may be added, but conflicts will be resolved
        in favor of the existing graph."""
        idgen = idGenerator if idGenerator is not None else idFactory("%s_U" % self.getName())

        # establish correct strand
        self_strand = coerce_enum(self.strand, Strand, field="strand")
        other_strand = coerce_enum(other.strand, Strand, field="strand")
        strand = self_strand
        if self_strand.value in VALID_STRANDS:
            if other_strand.value in VALID_STRANDS and self_strand != other_strand:
                raise ValueError("Cannot merge graphs with conflicting strands.")
        elif other_strand.value in VALID_STRANDS:
            strand = other_strand

        # update graph name
        if keepName:
            newName = self.getName()
        else:
            names = set(self.getName().split(","))
            names.add(other.getName())
            newName = ",".join(names)

        if mergeEnds:
            updateRoot(other, self)
            updateLeaf(other, self)

        result = SpliceGraph(newName, self.chromosome, strand.value)
        allNodes = self.resolvedNodes() + other.resolvedNodes()
        for node in allNodes:
            newNode = result.addNode(next(idgen), node.minpos, node.maxpos)
            for c in node.children:
                cNode = result.addNode(next(idgen), c.minpos, c.maxpos)
                newNode.addChild(cNode)
            for p in node.parents:
                pNode = result.addNode(next(idgen), p.minpos, p.maxpos)
                pNode.addChild(newNode)
            for iso in node.isoformSet:
                newNode.addIsoform(iso)
            for k in node.attrs:
                if k == AS_KEY:
                    for f in node.altForms():
                        newNode.addAltForm(f)
                elif k != ISO_KEY:
                    newNode.attrs[k] = node.attrs[k]
        return result

    def unresolvedNodes(self):
        """Returns a list of all unresolved nodes in the graph."""
        return [x for x in self.nodeDict.values() if x.isUnresolved()]

    def upstreamOf(self, a, b):
        """Returns true if a is upstream of b; false otherwise."""
        return (a < b) if self.strand == "+" else (b < a)

    def validate(self, halt: bool = False) -> str | None:
        """Returns None if the splicegraph is valid; otherwise returns a reason it is invalid."""
        reason: str | None = None
        all_nodes = list(self.nodeDict.values())
        ## allNodes   = self.resolvedNodes()
        node_set = set(all_nodes)
        node_list = list(node_set)
        all_node_ids = [node.id for node in all_nodes]
        unique_ids = [node.id for node in node_set]
        roots = self.getRoots()
        leaves = self.getLeaves()
        if len(roots) == 0 or len(leaves) == 0:
            reason = f"Graph is missing roots or leaves ({len(roots)} roots, {len(leaves)} leaves)"

        for n in all_nodes:
            if reason:
                break
            # Detect duplicate nodes: only one will appear in set
            if n.id not in unique_ids:
                other = node_list[node_list.index(n)]
                reason = (
                    f"Duplicate node {n.id} ({n.minpos}-{n.maxpos}) matches "
                    f"{other.id} ({other.minpos}-{other.maxpos})"
                )

            # Detect invalid child/parent ids
            for o in n.children:
                if o.id not in all_node_ids:
                    reason = f"Node {n.id}: child {o.id} not in graph"
                elif o.id not in unique_ids:
                    twin = node_list[node_list.index(o)]
                    reason = (
                        f"Node {n.id}: child {o.id} ({o.minpos}-{o.maxpos}) "
                        f"is not unique ({twin.id} {twin.minpos}-{twin.maxpos})"
                    )

            for o in n.parents:
                if o.id not in all_node_ids:
                    reason = f"Node {n.id} parent {o.id} not in graph"
                elif o.id not in unique_ids:
                    twin = node_list[node_list.index(o)]
                    reason = (
                        f"Node {n.id} parent {o.id} ({o.minpos}-{o.maxpos}) "
                        f"is not unique ({twin.id} {twin.minpos}-{twin.maxpos})"
                    )

        if reason and halt:
            LOGGER.error(
                "illegal_graph_validation_failed",
                nodes=[str(node) for node in all_nodes],
                reason=reason,
            )
            raise ValueError(f"Illegal graph:\n{reason}")
        elif reason:
            return reason
        return None

    def writeGFF(self, fileRef: str | TextIO, haltOnError: bool = False) -> bool:
        """Writes a splice graph to a file in GFF format.  The file may be given
        either as an output stream or as a file path.
        """
        reason = self.validate()
        if reason:
            if haltOnError:
                raise ValueError(
                    f'Cannot write invalid splice graph {self.getName()} to file:\n"{reason}"\n'
                )
            else:
                LOGGER.warning(
                    "writing_invalid_splice_graph",
                    graph_name=self.getName(),
                    reason=reason,
                )

        roots = set(self.getRoots())
        other = set(self.nodeDict.values()) - roots

        def _write_graph(out_stream: TextIO) -> None:
            out_stream.write(f"{self.gffString()}\n")
            for p in roots:
                out_stream.write(f"{p.gffString()}\n")
            for o in other:
                out_stream.write(f"{o.gffString()}\n")

        if isinstance(fileRef, str):
            with open(fileRef, "w") as out_stream:
                _write_graph(out_stream)
        else:
            _write_graph(fileRef)
        return True


class SpliceGraphParser:
    """Class that parses a GFF file filled with splice graphs and provides an
    iterator over each graph in the file."""

    def __init__(self, fileRef: str | TextIO, *, verbose: bool = False) -> None:
        """Parses a GFF file filled with splice graphs and returns each graph in an iterator.
        The file may be given either as an input stream or as a file path."""
        self.verbose = verbose

        if isinstance(fileRef, str):
            self.instream = ez_open(fileRef)
        else:
            self.instream = fileRef

        if self.instream is None:
            raise ValueError("No input file stream given.")

        self.graphDict: dict[str, SpliceGraph] = {}
        self._graph_keys: list[str] = []
        self.loadFromFile()

    def __iter__(self) -> SpliceGraphParser:
        """Iterator implementation."""
        return self

    def __next__(self) -> SpliceGraph:
        """Iterator implementation."""
        if self.graphId >= len(self._graph_keys):
            raise StopIteration
        key = self._graph_keys[self.graphId]
        self.graphId += 1
        return self.graphDict[key]

    def __len__(self) -> int:
        """Returns the number of nodes in the graph."""
        return len(self.graphDict)

    def loadFromFile(self) -> None:
        """Loads all graphs stored in a GFF file."""
        lineNo = 0
        graph: SpliceGraph | None = None
        # aliases keep track of nodes with different ids but identical start/end positions
        alias: dict[str, str] = {}
        edges: set[tuple[str, str]] = set()
        indicator = ProgressIndicator(100000, verbose=self.verbose)
        try:
            for line in self.instream:
                indicator.update()
                lineNo += 1
                if line.startswith("#"):
                    continue
                s = line.strip()
                parts = s.split("\t")

                try:
                    recType = coerce_enum(parts[2].lower(), RecordType, field="record_type").value
                    start = int(parts[3])
                    end = int(parts[4])
                except IndexError:
                    raise ValueError(
                        "Illegal record in splice graph file at line %d:\n\t%s" % (lineNo, s)
                    )
                except ValueError:
                    raise ValueError(
                        "Illegal record type in splice graph file at line %d:\n\t%s" % (lineNo, s)
                    )

                if recType not in VALID_RECTYPES:
                    raise ValueError(
                        "Illegal record type in splice graph file at line %d:\n\t%s" % (lineNo, s)
                    )

                attrs = self._parse_attributes(parts[-1], lineNo)

                try:
                    node_id = attrs[ID_ATTR]
                except KeyError:
                    raise ValueError(
                        "GFF attribute field '%s' has no ID at line %d" % (parts[-1], lineNo)
                    )

                if recType in VALID_GENES:
                    # Add edges to previous graph
                    if graph is not None:
                        for parent_id, child_id in edges:
                            graph.addEdge(alias[parent_id], alias[child_id])
                    # Start new graph
                    graph = SpliceGraph(name=node_id, chromosome=parts[0], strand=parts[6])
                    graph.minpos = min(start, end)
                    graph.maxpos = max(start, end)
                    for k in attrs:
                        if k not in KNOWN_ATTRS:
                            graph.attrs[k] = attrs[k]
                    self.graphDict[node_id] = graph
                    edges = set()
                    alias = {}
                elif graph is None:
                    raise ValueError("Graph feature found before graph header at line %d" % lineNo)
                else:
                    node = graph.addNode(node_id, start, end)
                    alias[node_id] = node.id
                    for k in attrs:
                        if k == AS_KEY:
                            node.addFormsFromString(attrs[k])
                        elif k == ISO_KEY and attrs[k]:
                            node.addIsoformString(attrs[k])
                        elif k in [START_CODON_KEY, END_CODON_KEY] and attrs[k]:
                            node.attrs[k] = {int(x) for x in attrs[k].split(",")}
                        elif k not in KNOWN_ATTRS:
                            node.addAttribute(k, attrs[k])

                    if PARENT_ATTR in attrs:
                        parents = attrs[PARENT_ATTR].split(",")
                        for parent in parents:
                            edges.add((parent, node_id))

            if graph is not None:
                for parent_id, child_id in edges:
                    graph.addEdge(alias[parent_id], alias[child_id])
        finally:
            indicator.finish()

        # Initialize iterator counter
        self.graphId = 0
        self._graph_keys = list(self.graphDict.keys())

    @staticmethod
    def _parse_attributes(field: str, line_no: int) -> dict[str, str]:
        """Convert `a=b;c=d` into key/value pairs while preserving '=' in values."""
        attrs: dict[str, str] = {}
        for pair in field.split(";"):
            if not pair:
                continue
            key, sep, value = pair.partition("=")
            if not sep or not key:
                raise ValueError(
                    "Illegal attribute field '%s' at line %d in GFF file." % (field, line_no)
                )
            attrs[key] = value
        return attrs
