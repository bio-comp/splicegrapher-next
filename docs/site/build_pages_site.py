"""Build the GitHub Pages site for SGN tutorial notebooks."""

from __future__ import annotations

import re
import shutil
import subprocess
import sys
import textwrap
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class NotebookTarget:
    """Notebook rendering contract for the Pages site."""

    notebook_path: str
    output_html: str
    title: str
    summary: str


NOTEBOOK_TARGETS: tuple[NotebookTarget, ...] = (
    NotebookTarget(
        notebook_path="docs/notebooks/01_legacy_tutorial_foundation.ipynb",
        output_html="01_legacy_tutorial_foundation.html",
        title="Legacy Tutorial Foundation",
        summary=(
            "Core SGN tutorial flow grounded in the original SpliceGrapher "
            "workflow, updated for Python 3 execution."
        ),
    ),
    NotebookTarget(
        notebook_path="docs/notebooks/02_modern_dataset_track.ipynb",
        output_html="02_modern_dataset_track.html",
        title="Modern Dataset Track",
        summary=("Plant and human dataset curation checks using the modern SGN manifest contract."),
    ),
)

HEAD_CLOSE_RE = re.compile(r"</head>", flags=re.IGNORECASE)
BODY_TAG_RE = re.compile(r"<body[^>]*>", flags=re.IGNORECASE)


def inject_theme_shell(raw_html: str, css_href: str, script_href: str, home_href: str) -> str:
    """Inject shared theme assets and top-bar controls into a rendered page."""
    head_match = HEAD_CLOSE_RE.search(raw_html)
    if head_match is None:
        raise ValueError("Rendered page missing </head> tag")

    head_snippet = textwrap.dedent(
        f"""
        <meta name="color-scheme" content="light dark">
        <link rel="stylesheet" href="{css_href}">
        <script defer src="{script_href}"></script>
        """
    ).strip()

    themed_html = (
        raw_html[: head_match.start()] + head_snippet + "\n" + raw_html[head_match.start() :]
    )

    body_match = BODY_TAG_RE.search(themed_html)
    if body_match is None:
        raise ValueError("Rendered page missing <body> tag")

    top_bar = textwrap.dedent(
        f"""
        <header class="sgn-topbar">
          <div class="sgn-topbar-inner">
            <a class="sgn-home-link" href="{home_href}">splicegrapher-next tutorials</a>
            <button class="sgn-theme-toggle" type="button" aria-label="Toggle theme"></button>
          </div>
        </header>
        """
    ).strip()

    return themed_html[: body_match.end()] + "\n" + top_bar + themed_html[body_match.end() :]


def render_index_html(targets: tuple[NotebookTarget, ...]) -> str:
    """Render the landing page with links to themed notebook pages."""
    cards = []
    for target in targets:
        cards.append(
            textwrap.dedent(
                f"""
                <article class="sgn-card">
                  <h2>{target.title}</h2>
                  <p>{target.summary}</p>
                  <a class="sgn-card-link" href="notebooks/{target.output_html}">Open notebook</a>
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


def _run_nbconvert(repo_root: Path, notebook_rel_path: str, output_html_name: str) -> None:
    notebook_path = repo_root / notebook_rel_path
    output_dir = repo_root / "site" / "notebooks"
    output_dir.mkdir(parents=True, exist_ok=True)

    command = [
        sys.executable,
        "-m",
        "jupyter",
        "nbconvert",
        "--to",
        "html",
        "--template",
        "lab",
        "--execute",
        notebook_rel_path,
        "--output",
        output_html_name,
        "--output-dir",
        "site/notebooks",
    ]
    subprocess.run(command, check=True, cwd=repo_root)

    if not notebook_path.exists():
        raise FileNotFoundError(f"Notebook not found: {notebook_path}")


def build_pages_site(repo_root: Path) -> None:
    """Build themed index + themed notebook pages under site/."""
    site_dir = repo_root / "site"
    site_dir.mkdir(parents=True, exist_ok=True)
    notebooks_dir = site_dir / "notebooks"
    notebooks_dir.mkdir(parents=True, exist_ok=True)

    assets_src = repo_root / "docs" / "site" / "assets"
    assets_dst = site_dir / "assets"
    assets_dst.mkdir(parents=True, exist_ok=True)

    for asset_name in ("pages.css", "theme.js"):
        shutil.copy2(assets_src / asset_name, assets_dst / asset_name)

    for target in NOTEBOOK_TARGETS:
        _run_nbconvert(repo_root, target.notebook_path, target.output_html)
        html_path = notebooks_dir / target.output_html
        themed_html = inject_theme_shell(
            raw_html=html_path.read_text(encoding="utf-8"),
            css_href="../assets/pages.css",
            script_href="../assets/theme.js",
            home_href="../index.html",
        )
        html_path.write_text(themed_html, encoding="utf-8")

    index_html = render_index_html(NOTEBOOK_TARGETS)
    (site_dir / "index.html").write_text(index_html, encoding="utf-8")


def main() -> None:
    repo_root = Path(__file__).resolve().parents[2]
    build_pages_site(repo_root=repo_root)


if __name__ == "__main__":
    main()
