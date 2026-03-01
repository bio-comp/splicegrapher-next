#!/usr/bin/env python3
"""Run reproducible benchmark comparisons for GFF loaders."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from SpliceGrapher.formats.polars_gff_benchmark import (
    DEFAULT_DATASET_SIZES,
    REQUIRED_REAL_DATASETS,
    evaluation_to_json_dict,
    evaluation_to_markdown,
    run_single_cycle_evaluation,
)


def _parse_sizes(raw_values: list[str] | None) -> dict[str, int]:
    if not raw_values:
        return dict(DEFAULT_DATASET_SIZES)

    result: dict[str, int] = {}
    for value in raw_values:
        if "=" not in value:
            raise ValueError(f"Invalid --size value '{value}'; expected label=count")

        label, raw_count = value.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError("Dataset label cannot be empty")

        count = int(raw_count)
        if count <= 0:
            raise ValueError("Dataset gene count must be >= 1")

        result[label] = count

    return result


def _parse_real_datasets(raw_values: list[str]) -> dict[str, Path]:
    result: dict[str, Path] = {}
    for value in raw_values:
        if "=" not in value:
            raise ValueError(
                f"Invalid --real-dataset value '{value}'; expected label=/path/to/file.gff3"
            )

        label, raw_path = value.split("=", 1)
        label = label.strip()
        if not label:
            raise ValueError("Real dataset label cannot be empty")

        path = Path(raw_path).expanduser()
        if not path.is_file():
            raise FileNotFoundError(f"Real dataset file not found: {path}")

        result[label] = path

    if len(result) < REQUIRED_REAL_DATASETS:
        raise ValueError(
            f"At least {REQUIRED_REAL_DATASETS} real datasets are required "
            "for final go/no-go decision"
        )

    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--synthetic-work-dir",
        type=Path,
        default=Path(".benchmarks") / "gff",
        help="Synthetic fixture output directory.",
    )
    parser.add_argument(
        "--iterations",
        type=int,
        default=3,
        help="Number of benchmark repetitions per loader/workload.",
    )
    parser.add_argument(
        "--exons-per-gene",
        type=int,
        default=3,
        help="Synthetic exon count per transcript.",
    )
    parser.add_argument(
        "--size",
        action="append",
        default=None,
        help="Synthetic dataset override as label=count; may be passed multiple times.",
    )
    parser.add_argument(
        "--real-dataset",
        action="append",
        default=[],
        help="Real dataset input as label=/absolute/or/relative/path.gff[3].",
    )
    parser.add_argument(
        "--skip-polars",
        action="store_true",
        help="Skip polars dataframe loader benchmark even if polars is installed.",
    )
    parser.add_argument(
        "--json-out",
        type=Path,
        default=Path("docs") / "testing" / "polars_gff_benchmark_results.json",
        help="JSON output file.",
    )
    parser.add_argument(
        "--markdown-out",
        type=Path,
        default=Path("docs") / "testing" / "polars_gff_benchmark_results.md",
        help="Markdown output file.",
    )
    args = parser.parse_args()

    sizes = _parse_sizes(args.size)
    real_datasets = _parse_real_datasets(args.real_dataset)

    evaluation = run_single_cycle_evaluation(
        synthetic_work_dir=args.synthetic_work_dir,
        real_datasets=real_datasets,
        synthetic_dataset_sizes=sizes,
        exons_per_gene=args.exons_per_gene,
        iterations=args.iterations,
        include_polars=not args.skip_polars,
    )

    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(
        json.dumps(evaluation_to_json_dict(evaluation), indent=2) + "\n",
        encoding="utf-8",
    )

    args.markdown_out.parent.mkdir(parents=True, exist_ok=True)
    args.markdown_out.write_text(
        evaluation_to_markdown(evaluation) + "\n",
        encoding="utf-8",
    )

    sys.stderr.write(f"Wrote JSON benchmark report: {args.json_out}\n")
    sys.stderr.write(f"Wrote markdown benchmark report: {args.markdown_out}\n")
    sys.stderr.write(
        f"Decision: {evaluation.decision.recommendation} ({evaluation.decision.decision})\n"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
