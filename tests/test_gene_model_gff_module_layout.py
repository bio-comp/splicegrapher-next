from __future__ import annotations

from SpliceGrapher.formats.parsers import gene_model_gff
from SpliceGrapher.formats.parsers import gene_model_gff_context as parser_context
from SpliceGrapher.formats.parsers import gene_model_gff_handlers as parser_handlers
from SpliceGrapher.formats.parsers import gene_model_gff_records as parser_records


def test_gene_model_gff_facade_preserves_loader_boundary() -> None:
    assert callable(gene_model_gff.load_gene_model_records)


def test_gene_model_gff_split_modules_expose_context_records_and_handlers() -> None:
    assert hasattr(parser_context, "ParseContext")
    assert hasattr(parser_context, "ParseStats")
    assert hasattr(parser_records, "ParsedRecord")
    assert hasattr(parser_records, "parse_record_line")
    assert hasattr(parser_handlers, "RECORD_HANDLERS")
    assert hasattr(parser_handlers, "resolve_parent")
