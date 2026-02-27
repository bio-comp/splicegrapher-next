"""Unit tests for Pages site builder helpers."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path
from types import ModuleType

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "docs" / "site" / "build_pages_site.py"


def _load_builder_module() -> ModuleType:
    spec = importlib.util.spec_from_file_location("build_pages_site", SCRIPT_PATH)
    if spec is None or spec.loader is None:
        raise RuntimeError("Unable to load docs/site/build_pages_site.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def test_notebook_html_name_uses_notebook_stem() -> None:
    module = _load_builder_module()
    notebook_path = "docs/notebooks/02_modern_dataset_track.ipynb"

    assert module.notebook_html_name(notebook_path) == "02_modern_dataset_track.html"


def test_render_index_html_contains_jupyter_book_links() -> None:
    module = _load_builder_module()

    page_html = module.render_index_html(module.NOTEBOOK_TARGETS)

    assert "<title>splicegrapher-next Tutorials</title>" in page_html
    assert "assets/pages.css" in page_html
    assert "assets/theme.js" in page_html
    assert "book/notebooks/01_legacy_tutorial_foundation.html" in page_html
    assert "book/notebooks/02_modern_dataset_track.html" in page_html


def test_render_book_toc_yaml_contains_expected_structure() -> None:
    module = _load_builder_module()

    toc_yaml = module.render_book_toc_yaml(module.NOTEBOOK_TARGETS)

    assert "format: jb-book" in toc_yaml
    assert "root: notebooks/01_legacy_tutorial_foundation" in toc_yaml
    assert "file: notebooks/02_modern_dataset_track" in toc_yaml
    assert "title: Modern Dataset Track" in toc_yaml


def test_render_book_config_yaml_uses_execute_off() -> None:
    module = _load_builder_module()

    config_yaml = module.render_book_config_yaml()

    assert "execute_notebooks: off" in config_yaml
    assert 'copyright: "2026"' in config_yaml
