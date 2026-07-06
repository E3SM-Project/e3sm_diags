"""Compare manual complete-run netCDF outputs against a baseline tree.

This module provides the manual command-line entrypoint for checking a dev
complete-run output directory against a known baseline directory. It exists
to preserve the current netCDF comparison workflow while moving the logic
out of pytest-style module scope and into a reusable, explicit, developer-run
tool that can summarize differences and optionally write PNG debug artifacts.

Usage
-----
Run this module with ``python -m tests.complete_run.compare``.

Use ``python -m tests.complete_run.compare --help`` for more details on the
available flags and their usage.

If the baselines need to be updated, run the complete-run workflow using
``python -m tests.complete_run.run`` on the main branch and update the
``DEFAULT_BASELINE_DIR`` constant in this module to point to the new baseline
directory. The new baseline directory should be committed to the repository
so that it is available to all developers and CI workflows.
"""

from __future__ import annotations

import argparse
import subprocess
from datetime import datetime
from pathlib import Path
from typing import Sequence

from e3sm_diags.logger import _setup_child_logger, _setup_root_logger
from tests.complete_run.helpers import (
    ComparisonIssue,
    ComparisonSummary,
    compare_netcdf_trees,
)
from tests.complete_run.params import DEFAULT_RESULTS_DIR

logger = _setup_child_logger(__name__)

# The default absolute and relative tolerances are set to match the current
# manual complete-run comparison workflow. These values can be overridden at
# runtime with the --atol and --rtol flags. NOTE: the default absolute tolerance
# is set to 0.0 because absolute tolerance is not used for
# floating-point comparison as it is highly sensitive to difference
# calculations (e.g., test-ref).
DEFAULT_ATOL = 0.0
DEFAULT_RTOL = 1e-5

# The default baseline directory is set to the current public E3SM Diags baseline
# TODO: Need to generate new baselines on main and update default path below.
DEFAULT_BASELINE_DIR = f"{DEFAULT_RESULTS_DIR}/<REPLACE-TIMESTAMP>-main-<REPLACE-HASH>"


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for manual netCDF comparison."""
    _setup_root_logger()

    parser = _build_parser()
    args = parser.parse_args(argv)
    compare_files, compare_values = _normalize_modes(args.mode)

    if not compare_files and not compare_values:
        raise ValueError("At least one compare mode must be selected.")

    _validate_compare_dirs(args.dev_dir, args.baseline_dir)

    diff_artifact_dir = None
    if args.write_diff_pngs:
        diff_artifact_dir = args.diff_artifact_dir or str(
            Path(args.dev_dir) / "compare-diffs"
        )

    summary = compare_netcdf_trees(
        dev_root=args.dev_dir,
        baseline_root=args.baseline_dir,
        atol=args.atol,
        rtol=args.rtol,
        compare_values=compare_values,
        diff_artifact_dir=diff_artifact_dir,
    )
    if not compare_files:
        summary.missing_dev_files.clear()
        summary.missing_baseline_files.clear()

    _render_summary(
        dev_dir=args.dev_dir, baseline_dir=args.baseline_dir, summary=summary
    )

    show = set(args.show or [])
    if "all" in show or "missing-files" in show:
        _render_issue_details("Missing dev files", summary.missing_dev_files)
        _render_issue_details("Missing baseline files", summary.missing_baseline_files)
    if "all" in show or "missing-vars" in show:
        _render_issue_details("Missing variables", summary.missing_variables)
    if "all" in show or "nan-mismatches" in show:
        _render_issue_details(
            "NaN-location mismatches", summary.nan_location_mismatches
        )
    if "all" in show or "shape-mismatches" in show:
        _render_issue_details("Shape mismatches", summary.shape_mismatches)
    if "all" in show or "tolerance-failures" in show:
        _render_issue_details("Tolerance failures", summary.tolerance_failures)

    return 0 if summary.failure_count == 0 else 1


def _build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for netCDF comparison."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compare manual complete-run netCDF outputs against a baseline tree."
        ),
        epilog=(
            "Example: python -m tests.complete_run.compare "
            "--dev-dir </path/to/dev>"
            f"--baseline-dir {DEFAULT_BASELINE_DIR} "
            "--show missing-files --show tolerance-failures"
        ),
    )
    parser.add_argument(
        "--dev-dir",
        help="The dev directory to compare (required).",
    )
    parser.add_argument(
        "--baseline-dir",
        default=DEFAULT_BASELINE_DIR,
        help=f"The baseline directory to compare against (default: {DEFAULT_BASELINE_DIR}).",
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=DEFAULT_ATOL,
        help=f"Absolute tolerance for netCDF value comparison (default: {DEFAULT_ATOL}). ",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=DEFAULT_RTOL,
        help=f"Relative tolerance for netCDF value comparison (default: {DEFAULT_RTOL}).",
    )
    parser.add_argument(
        "--mode",
        action="append",
        default=["all"],
        choices=["all", "files", "data"],
        help=(
            "Comparison mode. Use files for tree matching only, data for shared "
            "netCDF value comparison, or all for both (default: all)."
        ),
    )
    parser.add_argument(
        "--show",
        action="append",
        default=["all"],
        choices=[
            "all",
            "missing-files",
            "missing-vars",
            "nan-mismatches",
            "shape-mismatches",
            "tolerance-failures",
        ],
        help="Optional detail sections to emit after the top-level summary (default: all).",
    )
    parser.add_argument(
        "--write-diff-pngs",
        action="store_true",
        default=False,
        help="Write PNG diff artifacts for mismatched shared files (default: False).",
    )
    parser.add_argument(
        "--diff-artifact-dir",
        default=None,
        help="Directory for optional diff PNGs. Defaults to <dev-dir>/compare-diffs (default: None).",
    )

    return parser


def _default_dev_dir() -> str:
    """Build the default dev comparison path used by the current manual flow."""
    timestamp = datetime.now().strftime("%y-%m-%d")
    branch_name = _get_git_branch_name()

    return (
        f"/global/cfs/cdirs/e3sm/www/e3sm_diags/complete_run/{timestamp}-{branch_name}"
    )


def _get_git_branch_name() -> str:
    """Get the current git branch name."""
    try:
        return (
            subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                stderr=subprocess.DEVNULL,
            )
            .strip()
            .decode("utf-8")
        )
    except subprocess.CalledProcessError:
        return "unknown"


def _normalize_modes(modes: list[str] | None) -> tuple[bool, bool]:
    """Normalize CLI compare modes into file and value comparison toggles."""
    if not modes or "all" in modes:
        return True, True

    compare_files = "files" in modes
    compare_values = "data" in modes

    return compare_files, compare_values


def _validate_compare_dirs(dev_dir: str, baseline_dir: str) -> None:
    """Validate the directory inputs for the compare workflow.

    Raises
    ------
    FileNotFoundError
        If either comparison root directory does not exist.
    """
    missing_dirs = [
        path_str for path_str in (dev_dir, baseline_dir) if not Path(path_str).exists()
    ]

    if missing_dirs:
        missing_message = "\n".join(missing_dirs)
        raise FileNotFoundError(
            f"One or more compare directories do not exist:\n{missing_message}"
        )


def _render_summary(
    *,
    dev_dir: str,
    baseline_dir: str,
    summary: ComparisonSummary,
) -> None:
    """Emit the concise top-level comparison summary."""
    logger.info("Complete-run comparison summary")
    logger.info("dev: %s", dev_dir)
    logger.info("baseline: %s", baseline_dir)
    logger.info("matched files: %s", len(summary.matching_files))
    logger.info("missing dev files: %s", len(summary.missing_dev_files))
    logger.info("missing baseline files: %s", len(summary.missing_baseline_files))
    logger.info("missing vars: %s", len(summary.missing_variables))
    logger.info("NaN-location mismatches: %s", len(summary.nan_location_mismatches))
    logger.info("shape mismatches: %s", len(summary.shape_mismatches))
    logger.info("tolerance failures: %s", len(summary.tolerance_failures))
    logger.info("shared files compared: %s", summary.compared_file_count)


def _render_issue_details(
    title: str,
    issues: Sequence[Path | ComparisonIssue],
) -> None:
    """Emit a detailed list of paths or structured comparison issues."""
    if not issues:
        return

    logger.info("%s", title)
    for issue in issues:
        if isinstance(issue, Path):
            logger.info("  %s", issue)
            continue

        message = f"  {issue.relative_path}"
        if issue.var_key is not None:
            message += f" [var={issue.var_key}]"
        if issue.detail:
            message += f": {issue.detail}"
        if issue.artifact_path is not None:
            message += f" | diff_png={issue.artifact_path}"
        logger.info("%s", message)


if __name__ == "__main__":
    raise SystemExit(main())
