"""Typed TOML configuration loading backed by pydantic-settings."""

from __future__ import annotations

from pathlib import Path
from typing import cast

from pydantic import BaseModel, Field, ValidationError
from pydantic_settings import (
    BaseSettings,
    PydanticBaseSettingsSource,
    SettingsConfigDict,
    TomlConfigSettingsSource,
)

DEFAULT_CONFIG = Path("splicegrapher.toml")
DEFAULT_ENV_PREFIX = "SGN_"


class ConfigError(ValueError):
    """Base configuration error."""


class ConfigNotFoundError(ConfigError):
    """Raised when a configuration file path does not exist."""


class ConfigValidationError(ConfigError):
    """Raised when configuration content fails schema validation."""


class UnsupportedConfigFormatError(ConfigError):
    """Raised when a configuration file extension is not supported."""


class SpliceGrapherPaths(BaseModel):
    """Core path configuration for gene model and FASTA reference."""

    gene_model: Path | None = None
    fasta_reference: Path | None = None


class SpliceGrapherConfig(BaseSettings):
    """Canonical settings model for splicegrapher-next."""

    model_config = SettingsConfigDict(
        toml_file=str(DEFAULT_CONFIG),
        env_prefix=DEFAULT_ENV_PREFIX,
        env_nested_delimiter="__",
        extra="forbid",
        case_sensitive=False,
        validate_default=True,
    )

    splicegrapher: SpliceGrapherPaths = Field(...)

    @property
    def gene_model(self) -> Path | None:
        return self.splicegrapher.gene_model

    @property
    def fasta_reference(self) -> Path | None:
        return self.splicegrapher.fasta_reference

    @classmethod
    def settings_customise_sources(
        cls,
        settings_cls: type[BaseSettings],
        init_settings: PydanticBaseSettingsSource,
        env_settings: PydanticBaseSettingsSource,
        dotenv_settings: PydanticBaseSettingsSource,
        file_secret_settings: PydanticBaseSettingsSource,
    ) -> tuple[PydanticBaseSettingsSource, ...]:
        return (
            env_settings,
            TomlConfigSettingsSource(settings_cls),
            init_settings,
            dotenv_settings,
            file_secret_settings,
        )


def load_config(path: str | Path = DEFAULT_CONFIG) -> SpliceGrapherConfig:
    """Load canonical TOML config with env override support."""
    config_path = Path(path)
    if config_path.suffix.lower() != ".toml":
        raise UnsupportedConfigFormatError(
            f"Unsupported configuration format for {config_path!s}; expected .toml"
        )

    if not config_path.exists() or not config_path.is_file():
        raise ConfigNotFoundError(f"Configuration file not found at {config_path!s}")

    config_class = _config_class_for_toml_file(config_path)
    settings_loader = cast(type[BaseSettings], config_class)
    try:
        return cast(SpliceGrapherConfig, settings_loader())
    except ValidationError as exc:
        raise ConfigValidationError(
            f"Invalid TOML configuration values in {config_path!s}: {exc}"
        ) from exc


def _config_class_for_toml_file(config_path: Path) -> type[SpliceGrapherConfig]:
    model_config = SpliceGrapherConfig.model_config.copy()
    model_config["toml_file"] = str(config_path)
    return cast(
        type[SpliceGrapherConfig],
        type(
            "SpliceGrapherFileConfig",
            (SpliceGrapherConfig,),
            {"model_config": model_config},
        ),
    )


__all__ = [
    "ConfigError",
    "ConfigNotFoundError",
    "ConfigValidationError",
    "DEFAULT_CONFIG",
    "DEFAULT_ENV_PREFIX",
    "SpliceGrapherConfig",
    "SpliceGrapherPaths",
    "UnsupportedConfigFormatError",
    "load_config",
]
