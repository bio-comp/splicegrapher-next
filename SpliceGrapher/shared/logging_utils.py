"""Shared structlog bootstrap helpers."""

from __future__ import annotations

import sys
from typing import cast

import structlog
from structlog.typing import FilteringBoundLogger

_CONFIGURED = False


def configure_logging(*, level: int = 20) -> None:
    """Configure structlog once for SGN modules/scripts."""
    global _CONFIGURED
    if _CONFIGURED:
        return

    timestamper = structlog.processors.TimeStamper(fmt="iso", utc=True)
    structlog.configure(
        processors=[
            structlog.contextvars.merge_contextvars,
            structlog.processors.add_log_level,
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


def get_logger(
    name: str | None = None,
    *,
    level: int | None = None,
) -> FilteringBoundLogger:
    """Return a configured structlog logger."""
    if not _CONFIGURED:
        configure_logging(level=20 if level is None else level)
    return cast(FilteringBoundLogger, structlog.get_logger(name))
