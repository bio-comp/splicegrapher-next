from __future__ import annotations

from SpliceGrapher.core.enums import AlternativeSplicingEvent
from SpliceGrapher.core.splice_graph.graph import SpliceGraph
from SpliceGrapher.core.splicing_events import annotate_graph_events


def test_annotate_graph_events_marks_skipped_exon_for_bypass_motif() -> None:
    graph = SpliceGraph("skip", "chr1", "+")
    graph.add_node("a", 10, 20)
    graph.add_node("b", 30, 40)
    graph.add_node("c", 50, 60)
    graph.add_edge("a", "b")
    graph.add_edge("b", "c")
    graph.add_edge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ES in graph.node_dict["b"].alt_form_set


def test_annotate_graph_events_marks_retained_intron_for_containing_node() -> None:
    graph = SpliceGraph("retain", "chr1", "+")
    graph.add_node("a", 10, 20)
    graph.add_node("container", 15, 35)
    graph.add_node("c", 30, 40)
    graph.add_edge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.IR in graph.node_dict["container"].alt_form_set


def test_annotate_graph_events_marks_alt3_for_shared_parent_children() -> None:
    graph = SpliceGraph("alt3", "chr1", "+")
    graph.add_node("a", 10, 20)
    graph.add_node("b", 30, 40)
    graph.add_node("c", 35, 45)
    graph.add_edge("a", "b")
    graph.add_edge("a", "c")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALT3 in graph.node_dict["b"].alt_form_set
    assert AlternativeSplicingEvent.ALT3 in graph.node_dict["c"].alt_form_set


def test_annotate_graph_events_marks_alt5_for_shared_child_parents() -> None:
    graph = SpliceGraph("alt5", "chr1", "+")
    graph.add_node("b", 30, 40)
    graph.add_node("c", 35, 45)
    graph.add_node("d", 60, 70)
    graph.add_edge("b", "d")
    graph.add_edge("c", "d")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALT5 in graph.node_dict["b"].alt_form_set
    assert AlternativeSplicingEvent.ALT5 in graph.node_dict["c"].alt_form_set


def test_annotate_graph_events_upgrades_alt5_to_altb5_for_unique_acceptors() -> None:
    graph = SpliceGraph("altb5", "chr1", "+")
    graph.add_node("root", 1, 5)
    graph.add_node("alt_a", 10, 25)
    graph.add_node("alt_b", 15, 30)
    graph.add_node("acceptor_x", 40, 50)
    graph.add_node("acceptor_y", 60, 70)
    graph.add_edge("root", "alt_a")
    graph.add_edge("root", "alt_b")
    graph.add_edge("alt_a", "acceptor_x")
    graph.add_edge("alt_b", "acceptor_y")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALTB5 in graph.node_dict["alt_a"].alt_form_set
    assert AlternativeSplicingEvent.ALTB5 in graph.node_dict["alt_b"].alt_form_set
    assert AlternativeSplicingEvent.ALT5 not in graph.node_dict["alt_a"].alt_form_set
    assert AlternativeSplicingEvent.ALT5 not in graph.node_dict["alt_b"].alt_form_set


def test_annotate_graph_events_upgrades_alt3_to_altb3_for_unique_donors() -> None:
    graph = SpliceGraph("altb3", "chr1", "+")
    graph.add_node("donor_a", 1, 5)
    graph.add_node("donor_b", 2, 8)
    graph.add_node("alt_a", 10, 25)
    graph.add_node("alt_b", 15, 30)
    graph.add_edge("donor_a", "alt_a")
    graph.add_edge("donor_b", "alt_b")

    annotate_graph_events(graph)

    assert AlternativeSplicingEvent.ALTB3 in graph.node_dict["alt_a"].alt_form_set
    assert AlternativeSplicingEvent.ALTB3 in graph.node_dict["alt_b"].alt_form_set
    assert AlternativeSplicingEvent.ALT3 not in graph.node_dict["alt_a"].alt_form_set
    assert AlternativeSplicingEvent.ALT3 not in graph.node_dict["alt_b"].alt_form_set
