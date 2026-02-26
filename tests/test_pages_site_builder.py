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


def test_inject_theme_shell_adds_assets_and_navigation() -> None:
    module = _load_builder_module()
    raw_html = "<html><head><title>Notebook</title></head><body><main>Body</main></body></html>"

    themed = module.inject_theme_shell(
        raw_html=raw_html,
        css_href="../assets/pages.css",
        script_href="../assets/theme.js",
        home_href="../index.html",
    )

    assert 'href="../assets/pages.css"' in themed
    assert 'src="../assets/theme.js"' in themed
    assert 'class="sgn-topbar"' in themed
    assert 'href="../index.html"' in themed


def test_render_index_html_contains_notebook_links() -> None:
    module = _load_builder_module()

    page_html = module.render_index_html(module.NOTEBOOK_TARGETS)

    assert "<title>splicegrapher-next Tutorials</title>" in page_html
    assert "assets/pages.css" in page_html
    assert "assets/theme.js" in page_html
    assert "notebooks/01_legacy_tutorial_foundation.html" in page_html
    assert "notebooks/02_modern_dataset_track.html" in page_html
