"""Build the GitHub Pages site for SGN tutorial notebooks."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
import textwrap
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class NotebookTarget:
    """Notebook rendering contract for the Pages site."""

    notebook_path: str
    title: str
    summary: str


NOTEBOOK_TARGETS: tuple[NotebookTarget, ...] = (
    NotebookTarget(
        notebook_path="docs/notebooks/01_legacy_tutorial_foundation.ipynb",
        title="Paper-Lineage Tutorial Foundation",
        summary=(
            "Core SGN tutorial flow grounded in the published/original "
            "SpliceGrapher tutorial workflow, updated for Python 3 execution."
        ),
    ),
    NotebookTarget(
        notebook_path="docs/notebooks/02_modern_dataset_track.ipynb",
        title="Modern Dataset Track",
        summary=("Plant and human dataset curation checks using the modern SGN manifest contract."),
    ),
)


def notebook_html_name(notebook_path: str) -> str:
    """Resolve a notebook HTML output name from a notebook path."""
    return f"{Path(notebook_path).stem}.html"


def render_index_html(targets: tuple[NotebookTarget, ...]) -> str:
    """Render the landing page with links to Jupyter Book notebook pages."""
    cards = []
    for target in targets:
        notebook_html = notebook_html_name(target.notebook_path)
        cards.append(
            textwrap.dedent(
                f"""
                <article class="sgn-card">
                  <h2>{target.title}</h2>
                  <p>{target.summary}</p>
                  <a class="sgn-card-link" href="book/notebooks/{notebook_html}">Open notebook</a>
                </article>
                """
            ).strip()
        )

    cards_html = "\n".join(cards)
    return textwrap.dedent(
        f"""\
        <!doctype html>
        <html lang="en">
        <head>
          <meta charset="utf-8">
          <meta name="viewport" content="width=device-width, initial-scale=1">
          <meta name="color-scheme" content="light dark">
          <title>splicegrapher-next Tutorials</title>
          <link rel="stylesheet" href="assets/pages.css">
          <script defer src="assets/theme.js"></script>
        </head>
        <body>
          <header class="sgn-topbar">
            <div class="sgn-topbar-inner">
              <a class="sgn-home-link" href="https://github.com/bio-comp/splicegrapher-next">
                splicegrapher-next
              </a>
              <button class="sgn-theme-toggle" type="button" aria-label="Toggle theme"></button>
            </div>
          </header>
          <main class="sgn-shell">
            <section class="sgn-hero">
              <p class="sgn-eyebrow">Hands-on docs + executable notebooks</p>
              <h1>SpliceGrapher-Next Tutorial Hub</h1>
              <p>
                Modernized tutorials for splice-graph workflows, designed as
                user-facing documentation and integration smoke validation.
              </p>
            </section>
            <section class="sgn-grid">
              {cards_html}
            </section>
          </main>
        </body>
        </html>
        """
    )


def render_book_config_yaml() -> str:
    """Render temporary Jupyter Book config used for Pages build."""
    return textwrap.dedent(
        """\
        title: splicegrapher-next Tutorials
        author: splicegrapher-next maintainers
        copyright: "2026"
        exclude_patterns: [".DS_Store", "Thumbs.db"]

        execute:
          execute_notebooks: off

        html:
          title: splicegrapher-next Notebook Docs
        """
    )


def render_book_toc_yaml(targets: tuple[NotebookTarget, ...]) -> str:
    """Render temporary Jupyter Book TOC for notebook pages."""
    if not targets:
        raise ValueError("At least one notebook target is required")

    root_stem = Path(targets[0].notebook_path).stem
    lines = ["format: jb-book", f"root: notebooks/{root_stem}"]

    if len(targets) > 1:
        lines.append("chapters:")
        for target in targets[1:]:
            stem = Path(target.notebook_path).stem
            lines.append(f"  - file: notebooks/{stem}")
            lines.append(f"    title: {target.title}")

    return "\n".join(lines) + "\n"


def build_notebook_book(
    repo_root: Path,
    site_dir: Path,
    targets: tuple[NotebookTarget, ...],
) -> None:
    """Build notebook pages via Jupyter Book and place under site/book."""
    with tempfile.TemporaryDirectory(prefix="sgn-jupyter-book-") as tmp_dir_name:
        tmp_dir = Path(tmp_dir_name)
        book_root = tmp_dir / "book"
        notebooks_dir = book_root / "notebooks"
        notebooks_dir.mkdir(parents=True, exist_ok=True)

        for target in targets:
            source_path = repo_root / target.notebook_path
            if not source_path.exists():
                raise FileNotFoundError(f"Notebook target not found: {source_path}")
            shutil.copy2(source_path, notebooks_dir / source_path.name)

        (book_root / "_config.yml").write_text(render_book_config_yaml(), encoding="utf-8")
        (book_root / "_toc.yml").write_text(render_book_toc_yaml(targets), encoding="utf-8")

        jupyter_book_exe = shutil.which("jupyter-book")
        if jupyter_book_exe is None:
            raise FileNotFoundError("Unable to locate `jupyter-book` executable in PATH")

        command = [jupyter_book_exe, "build", str(book_root)]
        subprocess.run(command, check=True, cwd=repo_root)

        built_html_dir = book_root / "_build" / "html"
        if not built_html_dir.exists():
            raise FileNotFoundError(f"Jupyter Book output missing: {built_html_dir}")

        destination = site_dir / "book"
        if destination.exists():
            shutil.rmtree(destination)
        shutil.copytree(built_html_dir, destination)


def build_pages_site(repo_root: Path) -> None:
    """Build hybrid Pages output: custom landing + Jupyter Book notebooks."""
    site_dir = repo_root / "site"
    site_dir.mkdir(parents=True, exist_ok=True)

    assets_src = repo_root / "docs" / "site" / "assets"
    assets_dst = site_dir / "assets"
    assets_dst.mkdir(parents=True, exist_ok=True)

    for asset_name in ("pages.css", "theme.js"):
        shutil.copy2(assets_src / asset_name, assets_dst / asset_name)

    index_html = render_index_html(NOTEBOOK_TARGETS)
    (site_dir / "index.html").write_text(index_html, encoding="utf-8")
    build_notebook_book(repo_root=repo_root, site_dir=site_dir, targets=NOTEBOOK_TARGETS)


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    build_pages_site(repo_root=repo_root)


if __name__ == "__main__":
    main()
