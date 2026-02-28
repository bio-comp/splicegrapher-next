from SpliceGrapher.formats.annotation_io import load_gene_models


def loadGeneModels(path, **args):
    """Load gene models via the annotation I/O backend."""
    return load_gene_models(path, **args)
