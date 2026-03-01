"""Tests for progress reporting helpers."""

from __future__ import annotations

from SpliceGrapher.shared import progress as progress_module


class _FakeStream:
    def __init__(self, tty: bool) -> None:
        self._tty = tty
        self.writes: list[str] = []

    def isatty(self) -> bool:
        return self._tty

    def write(self, text: str) -> int:
        self.writes.append(text)
        return len(text)

    def flush(self) -> None:
        return None


class _FakeTqdm:
    def __init__(self, **kwargs) -> None:
        self.kwargs = kwargs
        self.total_updates = 0
        self.closed = False

    def update(self, amount: int = 1) -> None:
        self.total_updates += amount

    def close(self) -> None:
        self.closed = True


def test_progress_indicator_uses_tqdm_when_interactive(monkeypatch) -> None:
    fake_stream = _FakeStream(tty=True)
    created: list[_FakeTqdm] = []

    def make_bar(**kwargs) -> _FakeTqdm:
        bar = _FakeTqdm(**kwargs)
        created.append(bar)
        return bar

    monkeypatch.setattr(progress_module, "tqdm", make_bar)
    monkeypatch.setattr(progress_module.sys, "stderr", fake_stream)

    indicator = progress_module.ProgressIndicator(10, description="load", verbose=True)
    indicator.update()
    indicator.update()
    indicator.finish()

    assert len(created) == 1
    assert created[0].kwargs["disable"] is False
    assert created[0].kwargs["desc"] == "load"
    assert created[0].total_updates == 2
    assert created[0].closed is True


def test_progress_indicator_suppresses_output_when_not_interactive(monkeypatch) -> None:
    fake_stream = _FakeStream(tty=False)
    created: list[_FakeTqdm] = []

    def make_bar(**kwargs) -> _FakeTqdm:
        bar = _FakeTqdm(**kwargs)
        created.append(bar)
        return bar

    monkeypatch.setattr(progress_module, "tqdm", make_bar)
    monkeypatch.setattr(progress_module.sys, "stderr", fake_stream)

    indicator = progress_module.ProgressIndicator(10, description="load", verbose=True)
    indicator.update()
    indicator.update()
    indicator.finish()

    assert len(created) == 1
    assert created[0].kwargs["disable"] is True
    assert created[0].total_updates == 2
    assert fake_stream.writes == []


def test_progress_indicator_suppresses_output_when_verbose_is_false(monkeypatch) -> None:
    fake_stream = _FakeStream(tty=True)
    created: list[_FakeTqdm] = []

    def make_bar(**kwargs) -> _FakeTqdm:
        bar = _FakeTqdm(**kwargs)
        created.append(bar)
        return bar

    monkeypatch.setattr(progress_module, "tqdm", make_bar)
    monkeypatch.setattr(progress_module.sys, "stderr", fake_stream)

    indicator = progress_module.ProgressIndicator(10, description="load", verbose=False)
    indicator.update()
    indicator.finish()

    assert len(created) == 1
    assert created[0].kwargs["disable"] is True
    assert created[0].total_updates == 0
    assert fake_stream.writes == []
