# Configuration Guide

`splicegrapher-next` uses canonical TOML configuration validated through
`pydantic-settings`.

## Canonical Config

Create `splicegrapher.toml`:

```toml
[splicegrapher]
gene_model = "/absolute/path/to/genes.gtf"
fasta_reference = "/absolute/path/to/reference.fa"
```

## Python Loader API

```python
from SpliceGrapher.shared.config import load_config

config = load_config("splicegrapher.toml")

print(config.gene_model)
print(config.fasta_reference)
```

## Environment Override Contract (Nextflow/Docker Ready)

The settings model is configured with:

- `env_prefix = "SGN_"`
- `env_nested_delimiter = "__"`

This means orchestrators can override nested TOML fields without writing new
TOML files per run.

Example:

```bash
export SGN_SPLICEGRAPHER__GENE_MODEL=/mnt/data/genes.gtf
python -c "from SpliceGrapher.shared.config import load_config; print(load_config().gene_model)"
```

Source priority is:

1. Environment variables
2. TOML file
3. Code defaults (if defined)
