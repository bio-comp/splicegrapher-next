from __future__ import annotations

import importlib

from SpliceGrapher.formats.parsers import SpliceGraphParser, load_gene_model_records

parser_boundary = importlib.import_module("SpliceGrapher.formats.parsers.gene_model_gff")
parser_context = importlib.import_module("SpliceGrapher.formats.parsers.gene_model_gff_context")
parser_handlers = importlib.import_module("SpliceGrapher.formats.parsers.gene_model_gff_handlers")
parser_records = importlib.import_module("SpliceGrapher.formats.parsers.gene_model_gff_records")


def test_parsers_package_exposes_public_entrypoints() -> None:
    assert callable(load_gene_model_records)
    assert SpliceGraphParser.__name__ == "SpliceGraphParser"


def test_gene_model_gff_split_modules_expose_context_records_and_handlers() -> None:
    assert callable(parser_boundary.load_gene_model_records)
    assert hasattr(parser_context, "ParseContext")
    assert hasattr(parser_context, "ParseStats")
    assert hasattr(parser_records, "ParsedRecord")
    assert hasattr(parser_records, "parse_record_line")
    assert hasattr(parser_handlers, "RECORD_HANDLERS")
    assert hasattr(parser_handlers, "resolve_parent")
