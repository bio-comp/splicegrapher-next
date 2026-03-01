"""Behavior checks for shared process utility helpers."""

from __future__ import annotations

import subprocess
from typing import Any

from SpliceGrapher.shared import process_utils


def test_run_command_redirects_to_devnull_by_default(monkeypatch) -> None:
    captured: dict[str, Any] = {}

    def fake_call(command: str, **kwargs: Any) -> int:
        captured["command"] = command
        captured["stdout"] = kwargs.get("stdout")
        captured["stderr"] = kwargs.get("stderr")
        return 0

    monkeypatch.setattr(subprocess, "call", fake_call)
    process_utils.runCommand("echo hello", exitOnError=False)

    assert captured["command"] == "echo hello"
    assert captured["stdout"] is subprocess.DEVNULL
    assert captured["stderr"] is subprocess.DEVNULL


def test_run_command_preserves_explicit_stream_overrides(monkeypatch) -> None:
    captured: dict[str, Any] = {}

    def fake_call(command: str, **kwargs: Any) -> int:
        captured["command"] = command
        captured["stdout"] = kwargs.get("stdout")
        captured["stderr"] = kwargs.get("stderr")
        return 0

    monkeypatch.setattr(subprocess, "call", fake_call)

    process_utils.runCommand(
        "echo hello",
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        exitOnError=False,
    )

    assert captured["stdout"] is subprocess.PIPE
    assert captured["stderr"] is subprocess.PIPE
