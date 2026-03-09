from __future__ import annotations

from SpliceGrapher.formats import model_domain, model_index, models


def test_models_aggregator_reexports_domain_and_index_symbols() -> None:
    assert models.Gene is model_domain.Gene
    assert models.Transcript is model_domain.Transcript
    assert models.BaseFeature is model_domain.BaseFeature
    assert models.ChromosomeGeneIndex is model_index.ChromosomeGeneIndex
    assert models.Chromosome is model_index.Chromosome
    assert models.feature_sort_key is model_index.feature_sort_key
