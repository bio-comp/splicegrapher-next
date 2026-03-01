"""Behavior checks for shared structlog bootstrap helpers."""

from __future__ import annotations

from SpliceGrapher.shared import logging_utils


def test_configure_logging_is_idempotent() -> None:
    logging_utils._CONFIGURED = False
    logging_utils.configure_logging()
    logging_utils.configure_logging()

    assert logging_utils._CONFIGURED is True


def test_get_logger_returns_logger_with_info_method() -> None:
    logger = logging_utils.get_logger(__name__)
    assert hasattr(logger, "info")
