"""Focused tests for shared config/environment helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from SpliceGrapher.shared.config_utils import configMap, getEnvironmentValue


def test_config_map_reads_valid_ini_with_case_sensitive_keys(tmp_path: Path) -> None:
    cfg_path = tmp_path / "sample.ini"
    cfg_path.write_text("[SectionA]\nMixedCase = Value\nplain = x\n")

    parsed = configMap(cfg_path)

    assert parsed == {"SectionA": {"MixedCase": "Value", "plain": "x"}}


def test_config_map_raises_on_invalid_ini(tmp_path: Path) -> None:
    cfg_path = tmp_path / "broken.ini"
    cfg_path.write_text("[Broken\nmissing_close = 1\n")

    with pytest.raises(ValueError, match="Invalid configuration file"):
        configMap(cfg_path)


def test_config_map_raises_on_missing_file(tmp_path: Path) -> None:
    missing = tmp_path / "missing.ini"
    with pytest.raises(FileNotFoundError):
        configMap(missing)


def test_config_map_raises_on_invalid_option_value(tmp_path: Path) -> None:
    cfg_path = tmp_path / "invalid_value.ini"
    cfg_path.write_text("[SectionA]\nkey = %(missing)s\n")

    with pytest.raises(ValueError, match="Invalid configuration value"):
        configMap(cfg_path)


def test_get_environment_value_prefers_env_and_falls_back_to_default(monkeypatch) -> None:
    key = "SGN_CONFIG_UTILS_TEST"
    monkeypatch.delenv(key, raising=False)
    assert getEnvironmentValue(key, "fallback") == "fallback"

    monkeypatch.setenv(key, "present")
    assert getEnvironmentValue(key, "fallback") == "present"
