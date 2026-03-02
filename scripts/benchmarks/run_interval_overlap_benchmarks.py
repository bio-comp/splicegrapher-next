#!/usr/bin/env python3
"""Run reproducible overlap parity/performance benchmarks."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from pathlib import Path
from time import perf_counter

from SpliceGrapher.core.interval_helpers import InMemoryIntervalIndex

DEFAULT_INTERVAL_COUNT = 2000
DEFAULT_QUERY_COUNT = 2000
DEFAULT_ITERATIONS = 3
DEFAULT_MAX_SLOWDOWN_RATIO = 1.35


@dataclass(frozen=True, slots=True, order=True)
class _Interval:
    minpos: int
    maxpos: int


@dataclass(frozen=True, slots=True)
class OverlapBenchmarkResult:
    """Serializable overlap benchmark result payload."""

    interval_count: int
    query_count: int
    iterations: int
    legacy_mean_seconds: float
    indexed_mean_seconds: float
    indexed_slowdown_ratio: float
    parity_ok: bool
    max_slowdown_ratio: float

    @property
    def threshold_ok(self) -> bool:
        return self.parity_ok and self.indexed_slowdown_ratio <= self.max_slowdown_ratio


def _build_intervals(count: int) -> list[_Interval]:
    intervals: list[_Interval] = []
    for idx in range(count):
        start = (idx * 7) + 1
        end = start + 10 + (idx % 5)
        intervals.append(_Interval(minpos=start, maxpos=end))
    return intervals


def _build_queries(count: int) -> list[_Interval]:
    queries: list[_Interval] = []
    for idx in range(count):
        start = (idx * 5) + 3
        end = start + 8 + (idx % 3)
        queries.append(_Interval(minpos=start, maxpos=end))
    return queries


def _legacy_overlaps(intervals: list[_Interval], query: _Interval) -> list[_Interval]:
    return [
        interval
        for interval in intervals
        if interval.minpos <= query.maxpos and query.minpos <= interval.maxpos
    ]


def _run_legacy(
    intervals: list[_Interval],
    queries: list[_Interval],
    *,
    iterations: int,
) -> float:
    timings: list[float] = []
    for _ in range(iterations):
        start = perf_counter()
        checksum = 0
        for query in queries:
            checksum += len(_legacy_overlaps(intervals, query))
        _ = checksum
        timings.append(perf_counter() - start)
    return sum(timings) / len(timings)


def _run_indexed(
    intervals: list[_Interval],
    queries: list[_Interval],
    *,
    iterations: int,
) -> float:
    index = InMemoryIntervalIndex(intervals)
    timings: list[float] = []
    for _ in range(iterations):
        start = perf_counter()
        checksum = 0
        for query in queries:
            checksum += len(index.overlaps(query))
        _ = checksum
        timings.append(perf_counter() - start)
    return sum(timings) / len(timings)


def _parity_check(intervals: list[_Interval], queries: list[_Interval]) -> bool:
    index = InMemoryIntervalIndex(intervals)
    for query in queries:
        legacy = _legacy_overlaps(intervals, query)
        indexed = index.overlaps(query)
        if legacy != indexed:
            return False
    return True


def _to_markdown(result: OverlapBenchmarkResult) -> str:
    status = "pass" if result.threshold_ok else "fail"
    return "\n".join(
        [
            "# Interval Overlap Benchmark",
            "",
            f"- Status: **{status}**",
            f"- Parity: `{result.parity_ok}`",
            f"- Legacy mean seconds: `{result.legacy_mean_seconds:.6f}`",
            f"- Indexed mean seconds: `{result.indexed_mean_seconds:.6f}`",
            f"- Indexed slowdown ratio: `{result.indexed_slowdown_ratio:.3f}`",
            f"- Max slowdown ratio: `{result.max_slowdown_ratio:.3f}`",
            f"- Interval count: `{result.interval_count}`",
            f"- Query count: `{result.query_count}`",
            f"- Iterations: `{result.iterations}`",
        ]
    )


def run_benchmark(
    *,
    interval_count: int,
    query_count: int,
    iterations: int,
    max_slowdown_ratio: float,
) -> OverlapBenchmarkResult:
    if interval_count <= 0:
        raise ValueError("interval_count must be >= 1")
    if query_count <= 0:
        raise ValueError("query_count must be >= 1")
    if iterations <= 0:
        raise ValueError("iterations must be >= 1")
    if max_slowdown_ratio <= 0.0:
        raise ValueError("max_slowdown_ratio must be > 0.0")

    intervals = _build_intervals(interval_count)
    queries = _build_queries(query_count)
    legacy_mean = _run_legacy(intervals, queries, iterations=iterations)
    indexed_mean = _run_indexed(intervals, queries, iterations=iterations)
    slowdown = indexed_mean / legacy_mean if legacy_mean > 0.0 else float("inf")
    parity_ok = _parity_check(intervals, queries)

    return OverlapBenchmarkResult(
        interval_count=interval_count,
        query_count=query_count,
        iterations=iterations,
        legacy_mean_seconds=legacy_mean,
        indexed_mean_seconds=indexed_mean,
        indexed_slowdown_ratio=slowdown,
        parity_ok=parity_ok,
        max_slowdown_ratio=max_slowdown_ratio,
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--interval-count", type=int, default=DEFAULT_INTERVAL_COUNT)
    parser.add_argument("--query-count", type=int, default=DEFAULT_QUERY_COUNT)
    parser.add_argument("--iterations", type=int, default=DEFAULT_ITERATIONS)
    parser.add_argument(
        "--max-slowdown-ratio",
        type=float,
        default=DEFAULT_MAX_SLOWDOWN_RATIO,
        help="Fail threshold for indexed/legacy runtime ratio when --assert-threshold is set.",
    )
    parser.add_argument(
        "--assert-threshold",
        action="store_true",
        help="Exit non-zero when parity fails or slowdown ratio exceeds threshold.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=Path("docs") / "testing" / "interval_overlap_benchmark_results.json",
    )
    parser.add_argument(
        "--markdown-out",
        type=Path,
        default=Path("docs") / "testing" / "interval_overlap_benchmark_results.md",
    )
    args = parser.parse_args()

    result = run_benchmark(
        interval_count=args.interval_count,
        query_count=args.query_count,
        iterations=args.iterations,
        max_slowdown_ratio=args.max_slowdown_ratio,
    )

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(asdict(result), indent=2) + "\n", encoding="utf-8")

    args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_out.write_text(_to_markdown(result) + "\n", encoding="utf-8")

    sys.stderr.write(f"Wrote JSON benchmark report: {args.json_out}\n")
    sys.stderr.write(f"Wrote markdown benchmark report: {args.markdown_out}\n")
    sys.stderr.write(
        "Threshold check: "
        f"{'pass' if result.threshold_ok else 'fail'} "
        f"(ratio={result.indexed_slowdown_ratio:.3f}, limit={result.max_slowdown_ratio:.3f})\n"
    )

    if args.assert_threshold and not result.threshold_ok:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
