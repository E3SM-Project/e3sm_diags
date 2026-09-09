"""Run and validate a complete diagnostics candidate against ``latest-main``.

This module combines the expensive complete-run workflow with the repeatable
comparison step. It writes the immutable candidate results and comparison
artifacts before returning the comparison status, so expected changes remain
available for review.

Usage
-----
Run ``python -m tests.complete_run.validate`` on a compute node. Use
``--help`` for the complete-run and comparison options.
"""

from __future__ import annotations

import argparse
from typing import Sequence

from e3sm_diags.logger import _setup_child_logger, _setup_root_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from tests.complete_run import compare, run

logger = _setup_child_logger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    """Run a complete candidate and compare it with the accepted baseline."""
    _setup_root_logger()

    args = _build_parser().parse_args(argv)
    results = run._run_complete_run(args)
    _log_run_result(results)

    logger.info("Comparing complete-run results: %s", args.results_dir)
    return compare.main(_build_compare_argv(args))


def _build_parser() -> argparse.ArgumentParser:
    """Build a parser containing complete-run and comparison options."""
    parser = run._build_parser()
    parser.description = (
        "Run a complete diagnostics candidate and compare it with the accepted "
        "complete-run baseline."
    )
    parser.epilog = (
        "Example: python -m tests.complete_run.validate --set lat_lon --rtol 1e-5"
    )
    comparison_group = parser.add_argument_group("comparison options")
    comparison_group.add_argument(
        "--baseline-dir",
        default=compare.DEFAULT_BASELINE_DIR,
        help=(
            "Baseline directory to compare against "
            f"(default: {compare.DEFAULT_BASELINE_DIR})."
        ),
    )
    comparison_group.add_argument(
        "--atol",
        type=float,
        default=compare.DEFAULT_ATOL,
        help=f"Absolute tolerance for netCDF comparison (default: {compare.DEFAULT_ATOL}).",
    )
    comparison_group.add_argument(
        "--rtol",
        type=float,
        default=compare.DEFAULT_RTOL,
        help=f"Relative tolerance for netCDF comparison (default: {compare.DEFAULT_RTOL}).",
    )
    comparison_group.add_argument(
        "--image-mismatch-threshold",
        type=float,
        default=compare.DEFAULT_IMAGE_MISMATCH_THRESHOLD,
        help=(
            "Allowed mismatched-pixel fraction for PNG comparison "
            f"(default: {compare.DEFAULT_IMAGE_MISMATCH_THRESHOLD})."
        ),
    )
    comparison_group.add_argument(
        "--mode",
        action="append",
        choices=["all", "files", "data", "images"],
        help="Comparison mode; repeat to select multiple modes (default: all).",
    )
    comparison_group.add_argument(
        "--show",
        action="append",
        choices=[
            "all",
            "missing-files",
            "missing-vars",
            "nan-mismatches",
            "shape-mismatches",
            "tolerance-failures",
            "images",
        ],
        help="Comparison detail sections to log (default: all).",
    )
    comparison_group.add_argument(
        "--diff-artifact-dir",
        help="Directory for numerical and image diff artifacts.",
    )
    comparison_group.add_argument(
        "--report-dir",
        help="Directory for the comparison JSON report.",
    )
    return parser


def _build_compare_argv(args: argparse.Namespace) -> list[str]:
    """Convert combined CLI arguments to arguments accepted by ``compare``."""
    compare_argv = [
        "--dev-dir",
        args.results_dir,
        "--baseline-dir",
        str(args.baseline_dir),
        "--atol",
        str(args.atol),
        "--rtol",
        str(args.rtol),
        "--image-mismatch-threshold",
        str(args.image_mismatch_threshold),
        "--write-diff-pngs",
        "--write-diff-html",
    ]
    for mode in args.mode or ["all"]:
        compare_argv.extend(["--mode", mode])
    for section in args.show or ["all"]:
        compare_argv.extend(["--show", section])
    if args.diff_artifact_dir is not None:
        compare_argv.extend(["--diff-artifact-dir", args.diff_artifact_dir])
    if args.report_dir is not None:
        compare_argv.extend(["--report-dir", args.report_dir])

    return compare_argv


def _log_run_result(results: list[CoreParameter] | None) -> None:
    """Log the completed diagnostic result count without changing run semantics."""
    if results:
        logger.info("Complete-run finished with %s parameter results.", len(results))
    else:
        logger.warning(
            "Complete-run finished without explicit parameter results. "
            "Check the run log under the results directory if needed."
        )


if __name__ == "__main__":
    raise SystemExit(main())
