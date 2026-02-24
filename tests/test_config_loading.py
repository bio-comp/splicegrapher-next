from __future__ import annotations

from pathlib import Path

import pytest

from SpliceGrapher.shared.config import (
    ConfigNotFoundError,
    ConfigValidationError,
    SpliceGrapherConfig,
    UnsupportedConfigFormatError,
    load_config,
)


def test_load_config_parses_splicegrapher_section(tmp_path: Path) -> None:
    config_path = tmp_path / "splicegrapher.toml"
    config_path.write_text(
        "[splicegrapher]\ngene_model = '/tmp/genes.gtf'\nfasta_reference = '/tmp/reference.fa'\n",
        encoding="utf-8",
    )

    config = load_config(config_path)

    assert isinstance(config, SpliceGrapherConfig)
    assert config.gene_model == Path("/tmp/genes.gtf")
    assert config.fasta_reference == Path("/tmp/reference.fa")


def test_load_config_uses_nested_env_override(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
) -> None:
    config_path = tmp_path / "splicegrapher.toml"
    config_path.write_text("[splicegrapher]\ngene_model = '/tmp/default.gtf'\n", encoding="utf-8")
    monkeypatch.setenv("SGN_SPLICEGRAPHER__GENE_MODEL", "/tmp/from-env.gtf")

    config = load_config(config_path)

    assert config.gene_model == Path("/tmp/from-env.gtf")


def test_load_config_rejects_unsupported_extension(tmp_path: Path) -> None:
    bad_path = tmp_path / "splicegrapher.cfg"
    bad_path.write_text("[SpliceGrapher]\nGENE_MODEL = /tmp/legacy.gtf\n", encoding="utf-8")

    with pytest.raises(UnsupportedConfigFormatError):
        load_config(bad_path)


def test_load_config_requires_splicegrapher_section(tmp_path: Path) -> None:
    config_path = tmp_path / "splicegrapher.toml"
    config_path.write_text("project_name = 'wrong-schema'\n", encoding="utf-8")

    with pytest.raises(ConfigValidationError, match="splicegrapher"):
        load_config(config_path)


def test_load_config_validates_field_types(tmp_path: Path) -> None:
    config_path = tmp_path / "splicegrapher.toml"
    config_path.write_text("[splicegrapher]\ngene_model = 42\n", encoding="utf-8")

    with pytest.raises(ConfigValidationError, match="gene_model"):
        load_config(config_path)


def test_load_config_requires_existing_file(tmp_path: Path) -> None:
    with pytest.raises(ConfigNotFoundError):
        load_config(tmp_path / "missing.toml")
