"""Behavior checks for shared process utility helpers."""

from __future__ import annotations

import io
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


def test_log_message_emits_structured_event_and_mirrors_logstream(monkeypatch) -> None:
    captured: list[tuple[str, dict[str, Any]]] = []

    class StubLogger:
        def info(self, event: str, **kwargs: Any) -> None:
            captured.append((event, kwargs))

    monkeypatch.setattr(process_utils, "LOGGER", StubLogger())
    logstream = io.StringIO()

    process_utils.logMessage("hello world\n", logstream=logstream)

    assert captured == [("process_message", {"message": "hello world"})]
    assert logstream.getvalue() == "hello world\n"


def test_write_startup_message_emits_structured_event(monkeypatch) -> None:
    captured: list[tuple[str, dict[str, Any]]] = []

    class StubLogger:
        def info(self, event: str, **kwargs: Any) -> None:
            captured.append((event, kwargs))

    monkeypatch.setattr(process_utils, "LOGGER", StubLogger())
    monkeypatch.setattr(process_utils.sys, "argv", ["/tmp/example_script.py"])

    process_utils.writeStartupMessage()

    assert len(captured) == 1
    event, payload = captured[0]
    assert event == "startup"
    assert payload["script"] == "example_script.py"
