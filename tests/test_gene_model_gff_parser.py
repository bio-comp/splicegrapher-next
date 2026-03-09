from __future__ import annotations

import inspect

from SpliceGrapher.formats import gene_model as gm
from SpliceGrapher.formats.parsers import gene_model_gff as parser_boundary
from SpliceGrapher.formats.parsers import gene_model_gff_context as parser_context
from SpliceGrapher.formats.parsers import gene_model_gff_handlers as parser_handlers
from SpliceGrapher.formats.parsers import gene_model_gff_records as parser_records


def test_candidate_cache_uses_bounded_clear_policy() -> None:
    ctx = parser_context.ParseContext(
        model=gm.GeneModel(),
        require_notes=False,
        verbose=False,
        ignore_errors=False,
        parent_candidate_max_entries=2,
    )

    parser_handlers.candidate_parent_ids(ctx, "GENE1.1")
    parser_handlers.candidate_parent_ids(ctx, "GENE2.1")
    parser_handlers.candidate_parent_ids(ctx, "GENE3.1")

    assert len(ctx.parent_candidates) <= 2
    assert ctx.stats.parent_candidate_clear_count == 1


def test_parent_cache_uses_bounded_clear_policy() -> None:
    model = gm.GeneModel()
    model.add_chromosome(1, 100, "chr1")
    gene_one = gm.Gene("GENE1", None, 1, 10, "chr1", "+")
    gene_two = gm.Gene("GENE2", None, 20, 30, "chr1", "+")
    model.add_gene(gene_one)
    model.add_gene(gene_two)

    ctx = parser_context.ParseContext(
        model=model,
        require_notes=False,
        verbose=False,
        ignore_errors=False,
        parent_cache_max_entries=1,
    )

    assert parser_handlers.resolve_parent(ctx, "GENE1", "chr1") is gene_one
    assert parser_handlers.resolve_parent(ctx, "GENE2", "chr1") is gene_two

    assert len(ctx.parent_cache) == 1
    assert ctx.stats.parent_cache_clear_count == 1


def test_parser_modules_depend_on_models_module_not_gene_model_facade() -> None:
    source = "\n".join(
        [
            inspect.getsource(parser_boundary),
            inspect.getsource(parser_context),
            inspect.getsource(parser_records),
            inspect.getsource(parser_handlers),
        ]
    )

    assert "import SpliceGrapher.formats.models as model_domain" in source
    assert "import SpliceGrapher.formats.gene_model as gm" not in source
