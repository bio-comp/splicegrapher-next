"""Evaluate FASTA backend candidates without changing runtime package behavior."""

from __future__ import annotations

import argparse
import json
import random
import resource
import sys
import tempfile
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter
from typing import Callable

from pyfaidx import Fasta

from SpliceGrapher.formats import fasta as sg_fasta

try:  # pragma: no cover - optional dependency path
    import pyfastx
except ImportError:  # pragma: no cover - optional dependency path
    pyfastx = None


@dataclass(frozen=True, slots=True)
class Metric:
    mean_seconds: float
    max_seconds: float
    peak_mebibytes: float


@dataclass(frozen=True, slots=True)
class BackendResult:
    records: int
    iterate: Metric
    random_access: Metric


@dataclass(frozen=True, slots=True)
class IndexProbe:
    elapsed_seconds: float
    index_files: list[str]
    index_bytes: int


@dataclass(frozen=True, slots=True)
class EvaluationResult:
    config: dict[str, int]
    parity: dict[str, bool]
    contracts: dict[str, bool]
    sg: BackendResult
    pyfaidx: BackendResult
    pyfastx: BackendResult | None
    index_probe: dict[str, IndexProbe | None]
    recommendation: str


def _peak_rss_mebibytes() -> float:
    peak_rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    if sys.platform == "darwin":
        return peak_rss / (1024.0 * 1024.0)
    return peak_rss / 1024.0


def _write_fixture(path: Path, *, records: int, seq_len: int) -> list[str]:
    ids: list[str] = []
    with path.open("w", encoding="utf-8") as handle:
        for idx in range(records):
            record_id = f"seq{idx:04d} transcript {idx}"
            seq = ("ACGT" * (seq_len // 4 + 1))[:seq_len]
            handle.write(f">{record_id}\n{seq}\n")
            ids.append(record_id)
    return ids


def _measure(fn: Callable[[], object], *, iterations: int) -> Metric:
    timings: list[float] = []
    peaks: list[float] = []
    for _ in range(iterations):
        start = perf_counter()
        fn()
        timings.append(perf_counter() - start)
        peaks.append(_peak_rss_mebibytes())
    return Metric(
        mean_seconds=sum(timings) / len(timings),
        max_seconds=max(timings),
        peak_mebibytes=max(peaks),
    )


def _sg_iter(path: Path) -> dict[str, str]:
    return {rec.header: rec.sequence for rec in sg_fasta.FastaIterator(path)}


def _pyfaidx_iter(path: Path) -> dict[str, str]:
    fasta_file = Fasta(
        str(path),
        as_raw=True,
        read_long_names=True,
        sequence_always_upper=False,
    )
    try:
        return {key: str(fasta_file[key]) for key in fasta_file.keys()}
    finally:
        fasta_file.close()


def _pyfastx_iter(path: Path) -> dict[str, str]:
    if pyfastx is None:
        raise RuntimeError("pyfastx unavailable; run with `uv run --with pyfastx`.")
    fasta_file = pyfastx.Fasta(str(path), build_index=True, full_name=True)
    return {key: str(fasta_file[key]) for key in fasta_file.keys()}


def _sg_random_access(path: Path, ids: list[str], *, queries: int) -> None:
    for query_id in random.sample(ids, min(queries, len(ids))):
        record = sg_fasta.get_sequence(path, query_id)
        if record is None:
            raise RuntimeError(f"Missing sequence {query_id}")


def _pyfaidx_random_access(path: Path, ids: list[str], *, queries: int) -> None:
    fasta_file = Fasta(
        str(path),
        as_raw=True,
        read_long_names=True,
        sequence_always_upper=False,
    )
    try:
        for query_id in random.sample(ids, min(queries, len(ids))):
            _ = str(fasta_file[query_id])
    finally:
        fasta_file.close()


def _pyfastx_random_access(path: Path, ids: list[str], *, queries: int) -> None:
    if pyfastx is None:
        raise RuntimeError("pyfastx unavailable; run with `uv run --with pyfastx`.")
    fasta_file = pyfastx.Fasta(str(path), build_index=True, full_name=True)
    for query_id in random.sample(ids, min(queries, len(ids))):
        _ = str(fasta_file[query_id])


def _index_probe_pyfaidx(path: Path) -> IndexProbe:
    index_path = path.with_suffix(path.suffix + ".fai")
    if index_path.exists():
        index_path.unlink()
    start = perf_counter()
    fasta_file = Fasta(str(path), as_raw=True, read_long_names=True, sequence_always_upper=False)
    try:
        _ = next(iter(fasta_file.keys()))
    finally:
        fasta_file.close()
    elapsed = perf_counter() - start
    index_files = [index_path.name] if index_path.exists() else []
    index_bytes = index_path.stat().st_size if index_path.exists() else 0
    return IndexProbe(elapsed_seconds=elapsed, index_files=index_files, index_bytes=index_bytes)


def _index_probe_pyfastx(path: Path) -> IndexProbe | None:
    if pyfastx is None:
        return None
    index_path = Path(f"{path}.fxi")
    if index_path.exists():
        index_path.unlink()
    start = perf_counter()
    fasta_file = pyfastx.Fasta(str(path), build_index=True, full_name=True)
    _ = next(iter(fasta_file.keys()))
    elapsed = perf_counter() - start
    index_files = [index_path.name] if index_path.exists() else []
    index_bytes = index_path.stat().st_size if index_path.exists() else 0
    return IndexProbe(elapsed_seconds=elapsed, index_files=index_files, index_bytes=index_bytes)


def _supports_file_handle(path: Path) -> bool:
    if pyfastx is None:
        return False
    with path.open("r", encoding="utf-8") as handle:
        try:
            pyfastx.Fasta(handle, build_index=True)
        except Exception:
            return False
    return True


def evaluate(
    *,
    records: int,
    seq_len: int,
    iterations: int,
    random_access_queries: int,
) -> EvaluationResult:
    with tempfile.TemporaryDirectory() as temp_dir:
        temp_path = Path(temp_dir)
        fixture = temp_path / "probe.fa"
        ids = _write_fixture(fixture, records=records, seq_len=seq_len)

        sg_records = _sg_iter(fixture)
        pyfaidx_records = _pyfaidx_iter(fixture)
        parity: dict[str, bool] = {"pyfaidx_matches_sg": sg_records == pyfaidx_records}
        if pyfastx is not None:
            parity["pyfastx_matches_sg"] = sg_records == _pyfastx_iter(fixture)
        else:
            parity["pyfastx_matches_sg"] = False

        sg_iter_metric = _measure(lambda: _sg_iter(fixture), iterations=iterations)
        sg_access_metric = _measure(
            lambda: _sg_random_access(fixture, ids, queries=random_access_queries),
            iterations=iterations,
        )

        pyfaidx_iter_metric = _measure(lambda: _pyfaidx_iter(fixture), iterations=iterations)
        pyfaidx_access_metric = _measure(
            lambda: _pyfaidx_random_access(fixture, ids, queries=random_access_queries),
            iterations=iterations,
        )

        pyfastx_result: BackendResult | None = None
        if pyfastx is not None:
            pyfastx_iter_metric = _measure(lambda: _pyfastx_iter(fixture), iterations=iterations)
            pyfastx_access_metric = _measure(
                lambda: _pyfastx_random_access(fixture, ids, queries=random_access_queries),
                iterations=iterations,
            )
            pyfastx_result = BackendResult(
                records=len(sg_records),
                iterate=pyfastx_iter_metric,
                random_access=pyfastx_access_metric,
            )

        contracts = {
            "sg_supports_file_handles": True,
            "pyfastx_supports_file_handles": _supports_file_handle(fixture),
        }

        recommendation = "reject_default_keep_pyfaidx"
        if pyfastx is None:
            recommendation = "reject_default_pyfastx_not_installed"
        elif contracts["pyfastx_supports_file_handles"]:
            recommendation = "optional_backend_only"

        return EvaluationResult(
            config={
                "records": records,
                "seq_len": seq_len,
                "iterations": iterations,
                "random_access_queries": random_access_queries,
            },
            parity=parity,
            contracts=contracts,
            sg=BackendResult(
                records=len(sg_records),
                iterate=sg_iter_metric,
                random_access=sg_access_metric,
            ),
            pyfaidx=BackendResult(
                records=len(sg_records),
                iterate=pyfaidx_iter_metric,
                random_access=pyfaidx_access_metric,
            ),
            pyfastx=pyfastx_result,
            index_probe={
                "pyfaidx": _index_probe_pyfaidx(fixture),
                "pyfastx": _index_probe_pyfastx(fixture),
            },
            recommendation=recommendation,
        )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--records", type=int, default=3000)
    parser.add_argument("--seq-len", type=int, default=250)
    parser.add_argument("--iterations", type=int, default=5)
    parser.add_argument("--random-access-queries", type=int, default=500)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()

    result = evaluate(
        records=args.records,
        seq_len=args.seq_len,
        iterations=args.iterations,
        random_access_queries=args.random_access_queries,
    )
    payload = json.dumps(asdict(result), indent=2, sort_keys=True)
    if args.output is None:
        sys.stdout.write(f"{payload}\n")
        return
    args.output.write_text(payload, encoding="utf-8")
    sys.stdout.write(f"{args.output}\n")


if __name__ == "__main__":
    main()
