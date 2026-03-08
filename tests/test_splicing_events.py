from __future__ import annotations

from SpliceGrapher.core.enums import AlternativeSplicingEvent
from SpliceGrapher.core.splice_graph import SpliceGraph
from SpliceGrapher.core.splicing_events import annotate_graph_events


def test_annotate_graph_events_marks_skipped_exon_for_bypass_motif() -> None:
    graph = SpliceGraph("skip", "chr1", "+")
    graph.addNode("a", 10, 20)
    graph.addNode("b", 30, 40)
    graph.addNode("c", 50, 60)
    graph.addEdge("a", "b")
    graph.addEdge("b", "c")
    graph.addEdge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ES in graph.nodeDict["b"].altFormSet


def test_annotate_graph_events_marks_retained_intron_for_containing_node() -> None:
    graph = SpliceGraph("retain", "chr1", "+")
    graph.addNode("a", 10, 20)
    graph.addNode("container", 15, 35)
    graph.addNode("c", 30, 40)
    graph.addEdge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.IR in graph.nodeDict["container"].altFormSet


def test_annotate_graph_events_marks_alt3_for_shared_parent_children() -> None:
    graph = SpliceGraph("alt3", "chr1", "+")
    graph.addNode("a", 10, 20)
    graph.addNode("b", 30, 40)
    graph.addNode("c", 35, 45)
    graph.addEdge("a", "b")
    graph.addEdge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALT3 in graph.nodeDict["b"].altFormSet
    assert AlternativeSplicingEvent.ALT3 in graph.nodeDict["c"].altFormSet


def test_annotate_graph_events_marks_alt5_for_shared_child_parents() -> None:
    graph = SpliceGraph("alt5", "chr1", "+")
    graph.addNode("b", 30, 40)
    graph.addNode("c", 35, 45)
    graph.addNode("d", 60, 70)
    graph.addEdge("b", "d")
    graph.addEdge("c", "d")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALT5 in graph.nodeDict["b"].altFormSet
    assert AlternativeSplicingEvent.ALT5 in graph.nodeDict["c"].altFormSet


def test_annotate_graph_events_upgrades_alt5_to_altb5_for_unique_acceptors() -> None:
    graph = SpliceGraph("altb5", "chr1", "+")
    graph.addNode("root", 1, 5)
    graph.addNode("alt_a", 10, 25)
    graph.addNode("alt_b", 15, 30)
    graph.addNode("acceptor_x", 40, 50)
    graph.addNode("acceptor_y", 60, 70)
    graph.addEdge("root", "alt_a")
    graph.addEdge("root", "alt_b")
    graph.addEdge("alt_a", "acceptor_x")
    graph.addEdge("alt_b", "acceptor_y")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALTB5 in graph.nodeDict["alt_a"].altFormSet
    assert AlternativeSplicingEvent.ALTB5 in graph.nodeDict["alt_b"].altFormSet
    assert AlternativeSplicingEvent.ALT5 not in graph.nodeDict["alt_a"].altFormSet
    assert AlternativeSplicingEvent.ALT5 not in graph.nodeDict["alt_b"].altFormSet


def test_annotate_graph_events_upgrades_alt3_to_altb3_for_unique_donors() -> None:
    graph = SpliceGraph("altb3", "chr1", "+")
    graph.addNode("donor_a", 1, 5)
    graph.addNode("donor_b", 2, 8)
    graph.addNode("alt_a", 10, 25)
    graph.addNode("alt_b", 15, 30)
    graph.addEdge("donor_a", "alt_a")
    graph.addEdge("donor_b", "alt_b")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALTB3 in graph.nodeDict["alt_a"].altFormSet
    assert AlternativeSplicingEvent.ALTB3 in graph.nodeDict["alt_b"].altFormSet
    assert AlternativeSplicingEvent.ALT3 not in graph.nodeDict["alt_a"].altFormSet
    assert AlternativeSplicingEvent.ALT3 not in graph.nodeDict["alt_b"].altFormSet
