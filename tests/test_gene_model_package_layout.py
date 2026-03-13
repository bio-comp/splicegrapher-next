from __future__ import annotations

import SpliceGrapher.formats.gene_model.model as gene_model_model
import SpliceGrapher.formats.gene_model.repository as gene_model_repository
from SpliceGrapher.formats import gene_model


def test_gene_model_package_reexports_split_module_symbols() -> None:
    assert gene_model.GeneModel is gene_model_model.GeneModel
    assert gene_model.GeneModelRepository is gene_model_repository.GeneModelRepository
