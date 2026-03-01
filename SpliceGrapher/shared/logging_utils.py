"""Shared structlog bootstrap helpers."""

from __future__ import annotations

import sys
from typing import Any

import structlog

_CONFIGURED = False


def configure_logging(*, level: int = 20) -> None:
    """Configure structlog once for SGN modules/scripts."""
    global _CONFIGURED
    if _CONFIGURED:
        return

    timestamper = structlog.processors.TimeStamper(fmt="iso", utc=True)
    structlog.configure(
        processors=[
            structlog.stdlib.add_log_level,
            timestamper,
            structlog.processors.StackInfoRenderer(),
            structlog.dev.set_exc_info,
            structlog.dev.ConsoleRenderer(),
        ],
        wrapper_class=structlog.make_filtering_bound_logger(level),
        logger_factory=structlog.PrintLoggerFactory(file=sys.stderr),
        cache_logger_on_first_use=True,
    )
    _CONFIGURED = True


def get_logger(name: str | None = None) -> Any:
    """Return a configured structlog logger."""
    configure_logging()
    return structlog.get_logger(name)
