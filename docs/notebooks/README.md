# SGN Tutorial Notebooks

This directory contains executable tutorial notebooks for `splicegrapher-next`.

Current tracks:

- `01_legacy_tutorial_foundation.ipynb`
  - legacy tutorial framing with SGN modern manifest integration
- `02_modern_dataset_track.ipynb`
  - modern plant/human curated dataset contract checks

## Local Execution

```bash
uv sync --group dev
uv run jupyter nbconvert --to notebook --execute docs/notebooks/01_legacy_tutorial_foundation.ipynb --output /tmp/sgn-legacy.ipynb
uv run jupyter nbconvert --to notebook --execute docs/notebooks/02_modern_dataset_track.ipynb --output /tmp/sgn-modern.ipynb
```

JupyterLab interactive workflow:

```bash
uv run jupyter lab
```

## CI and Pages

- PR smoke execution is handled by `.github/workflows/notebook-smoke.yml`.
- GitHub Pages rendering/deploy is handled by
  `.github/workflows/docs-notebooks-pages.yml`.
