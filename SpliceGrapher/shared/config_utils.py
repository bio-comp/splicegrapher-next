"""Configuration and environment helpers extracted from shared.utils."""

from __future__ import annotations

import configparser
import os
from pathlib import Path


def configMap(cfgFile: str | Path) -> dict[str, dict[str, str]]:
    """Read an INI configuration file into a nested mapping.

    Keeps section and option names case-sensitive for compatibility with
    existing SGN config conventions.
    """
    cfg_path = cfgFile if isinstance(cfgFile, Path) else Path(cfgFile)
    if not cfg_path.is_file():
        raise FileNotFoundError(f"Configuration file not found: {cfg_path}")

    config = configparser.ConfigParser()
    config.optionxform = str
    try:
        with cfg_path.open("r", encoding="utf-8") as config_stream:
            config.read_file(config_stream)
    except (configparser.ParsingError, configparser.MissingSectionHeaderError) as exc:
        raise ValueError(f"Invalid configuration file {cfg_path}") from exc

    result: dict[str, dict[str, str]] = {}
    for sect in config.sections():
        values: dict[str, str] = {}
        for opt in config.options(sect):
            try:
                values[opt] = config.get(sect, opt)
            except configparser.Error as exc:
                raise ValueError(
                    f"Invalid configuration value for option '{opt}' in section '{sect}'"
                ) from exc
        result[sect] = values
    return result


def getEnvironmentValue(name: str, default: str | None = None) -> str | None:
    """Returns the value for the given environment variable name, if found;
    otherwise returns default value."""
    return os.getenv(name, default)
