from __future__ import annotations

import SpliceGrapher.formats.models.domain as model_domain
import SpliceGrapher.formats.models.index as model_index
from SpliceGrapher.formats import models


def test_models_package_reexports_domain_and_index_symbols() -> None:
    assert models.Gene is model_domain.Gene
    assert models.Transcript is model_domain.Transcript
    assert models.BaseFeature is model_domain.BaseFeature
    assert models.ChromosomeGeneIndex is model_index.ChromosomeGeneIndex
    assert models.Chromosome is model_index.Chromosome
    assert models.feature_sort_key is model_index.feature_sort_key
