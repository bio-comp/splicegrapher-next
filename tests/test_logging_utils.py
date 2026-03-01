"""Behavior checks for shared structlog bootstrap helpers."""

from __future__ import annotations

from typing import Any

import structlog

from SpliceGrapher.shared import logging_utils


def test_configure_logging_is_idempotent() -> None:
    logging_utils._CONFIGURED = False
    logging_utils.configure_logging()
    logging_utils.configure_logging()

    assert logging_utils._CONFIGURED is True


def test_get_logger_returns_logger_with_info_method() -> None:
    logger = logging_utils.get_logger(__name__)
    assert hasattr(logger, "info")


def test_configure_logging_uses_nonstdlib_processors(monkeypatch) -> None:
    captured: dict[str, Any] = {}

    def fake_configure(**kwargs: Any) -> None:
        captured.update(kwargs)

    logging_utils._CONFIGURED = False
    monkeypatch.setattr(structlog, "configure", fake_configure)

    logging_utils.configure_logging(level=10)

    processors = captured["processors"]
    assert processors[0] is structlog.contextvars.merge_contextvars
    assert processors[1] is structlog.processors.add_log_level


def test_get_logger_passes_level_to_first_configure_call(monkeypatch) -> None:
    seen: dict[str, Any] = {}
    sentinel_logger = structlog.get_logger("sentinel")

    def fake_configure_logging(*, level: int = 20) -> None:
        seen["level"] = level
        logging_utils._CONFIGURED = True

    logging_utils._CONFIGURED = False
    monkeypatch.setattr(logging_utils, "configure_logging", fake_configure_logging)
    monkeypatch.setattr(structlog, "get_logger", lambda name=None: sentinel_logger)

    logger = logging_utils.get_logger("sample", level=10)

    assert seen["level"] == 10
    assert logger is sentinel_logger
