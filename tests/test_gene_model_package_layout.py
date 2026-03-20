from __future__ import annotations

import SpliceGrapher.formats.gene_model.model as gene_model_model
from SpliceGrapher.formats import gene_model


def test_gene_model_package_reexports_only_runtime_symbols() -> None:
    assert gene_model.GeneModel is gene_model_model.GeneModel
    assert not hasattr(gene_model, "GeneModelRepository")
