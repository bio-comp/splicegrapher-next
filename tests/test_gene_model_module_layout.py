from __future__ import annotations

import SpliceGrapher.formats.models.constants as model_constants
import SpliceGrapher.formats.models.features as model_features
import SpliceGrapher.formats.models.gene as model_gene
import SpliceGrapher.formats.models.index as model_index
import SpliceGrapher.formats.models.locus as model_locus
import SpliceGrapher.formats.models.transcript as model_transcript
from SpliceGrapher.formats import models


def test_models_package_reexports_split_module_symbols() -> None:
    assert models.Gene is model_gene.Gene
    assert models.Transcript is model_transcript.Transcript
    assert models.BaseFeature is model_features.BaseFeature
    assert models.Locus is model_locus.Locus
    assert models.KNOWN_RECTYPES is model_constants.KNOWN_RECTYPES
    assert models.ChromosomeGeneIndex is model_index.ChromosomeGeneIndex
    assert models.Chromosome is model_index.Chromosome
    assert models.feature_sort_key is model_index.feature_sort_key
