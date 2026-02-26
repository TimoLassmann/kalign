"""Unified CLI entry point for downstream benchmarks.

Usage::

    python -m benchmarks.downstream --all -j 4
    python -m benchmarks.downstream --all -j 4 --quick
    python -m benchmarks.downstream calibration -j 4
    python -m benchmarks.downstream positive_selection --quick
    python -m benchmarks.downstream --figures -o benchmarks/figures/
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from pathlib import Path

logger = logging.getLogger("benchmarks.downstream")

PIPELINES = ["alignment_accuracy", "calibration", "positive_selection", "phylo_accuracy", "hmmer_detection"]


def _run_pipeline(name: str, params: dict) -> None:
    """Import and run a single pipeline by name."""
    if name == "alignment_accuracy":
        from .alignment_accuracy import run_pipeline
    elif name == "calibration":
        from .calibration import run_pipeline
    elif name == "positive_selection":
        from .positive_selection import run_pipeline
    elif name == "phylo_accuracy":
        from .phylo_accuracy import run_pipeline
    elif name == "hmmer_detection":
        from .hmmer_detection import run_pipeline
    else:
        raise ValueError(f"Unknown pipeline: {name!r}")

    logger.info("Starting pipeline: %s", name)
    start = time.perf_counter()
    run_pipeline(params)
    elapsed = time.perf_counter() - start
    logger.info("Pipeline %s completed in %.1f s", name, elapsed)


def _run_figures(results_dir: Path, output_dir: Path) -> None:
    """Generate all publication figures from saved results."""
    from .figures import generate_all_figures

    generate_all_figures(str(results_dir), str(output_dir))


def main(argv: list[str] | None = None) -> int:
    """Parse arguments and dispatch to pipelines."""
    parser = argparse.ArgumentParser(
        prog="python -m benchmarks.downstream",
        description="Kalign downstream application benchmarks",
    )

    parser.add_argument(
        "pipelines",
        nargs="*",
        default=[],
        metavar="PIPELINE",
        help="Pipeline(s) to run: alignment_accuracy, calibration, "
             "positive_selection, phylo_accuracy, hmmer_detection, or 'all'.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Run all four pipelines.",
    )
    parser.add_argument(
        "--figures",
        action="store_true",
        help="Generate publication figures from saved results.",
    )
    _default_jobs = os.cpu_count() or 1
    parser.add_argument(
        "-j", "--jobs",
        type=int,
        default=_default_jobs,
        help=f"Number of parallel jobs (default: {_default_jobs}, all cores). "
             "Use -j 1 for sequential execution with accurate timing.",
    )
    parser.add_argument(
        "--quick",
        action="store_true",
        help="Quick smoke test (5 cases per pipeline).",
    )
    parser.add_argument(
        "--data-dir",
        type=Path,
        default=Path("benchmarks/data"),
        help="Data directory (default: benchmarks/data).",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=Path("benchmarks/results"),
        help="Results directory (default: benchmarks/results).",
    )
    parser.add_argument(
        "-o", "--output-dir",
        type=Path,
        default=Path("benchmarks/figures"),
        help="Figure output directory (default: benchmarks/figures).",
    )
    parser.add_argument(
        "--methods",
        nargs="+",
        default=None,
        help="Override method list (default: pipeline-specific defaults).",
    )
    parser.add_argument(
        "--no-cache",
        action="store_true",
        help="Ignore cached results and recompute everything.",
    )
    parser.add_argument(
        "--full",
        action="store_true",
        help="Use full simulation grids and all methods (for publication).",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="Enable debug logging.",
    )
    parser.add_argument(
        "--depths",
        nargs="+",
        type=float,
        default=None,
        help="Filter simulations to specific tree depths (e.g., --depths 4.0).",
    )

    args = parser.parse_args(argv)

    # Configure logging
    level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(name)s %(levelname)s %(message)s",
        datefmt="%H:%M:%S",
    )

    # Validate pipeline names
    valid = set(PIPELINES + ["all"])
    for p in args.pipelines:
        if p not in valid:
            parser.error(f"Unknown pipeline: {p!r}. Choose from: {', '.join(PIPELINES)}, all")

    # Determine which pipelines to run
    pipelines_to_run: list[str] = []
    if args.all or "all" in args.pipelines:
        pipelines_to_run = list(PIPELINES)
    else:
        pipelines_to_run = list(args.pipelines)

    # Must specify at least one action
    if not pipelines_to_run and not args.figures:
        parser.print_help()
        return 1

    # Build common params dict
    params = {
        "data_dir": str(args.data_dir),
        "results_dir": str(args.results_dir),
        "n_jobs": args.jobs,
        "quick": args.quick,
        "no_cache": args.no_cache,
        "full": args.full,
    }
    if args.methods:
        params["methods"] = args.methods
    if args.depths:
        params["depths"] = args.depths

    # Run pipelines
    failed: list[str] = []
    for name in pipelines_to_run:
        try:
            _run_pipeline(name, params)
        except Exception:
            logger.exception("Pipeline %s failed", name)
            failed.append(name)

    # Generate figures
    if args.figures:
        try:
            _run_figures(args.results_dir, args.output_dir)
        except Exception:
            logger.exception("Figure generation failed")
            failed.append("figures")

    # Summary
    if failed:
        logger.error("Failed: %s", ", ".join(failed))
        return 1

    if pipelines_to_run:
        logger.info("All pipelines completed successfully.")
    if args.figures:
        logger.info("Figures saved to %s", args.output_dir)

    return 0


if __name__ == "__main__":
    sys.exit(main())
