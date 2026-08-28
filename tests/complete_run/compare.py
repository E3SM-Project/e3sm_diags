"""Compare manual complete-run netCDF and PNG outputs against a baseline tree.

This module provides the manual command-line entrypoint for checking a dev
complete-run output directory against a known baseline directory. It exists
to preserve the current complete-run comparison workflow while moving the logic
out of pytest-style module scope and into a reusable, explicit, developer-run
tool that can summarize numerical and visual differences and write PNG debug artifacts.
Each comparison also writes a machine-readable JSON report.

Usage
-----
Run this module with ``python -m tests.complete_run.compare``.

Use ``python -m tests.complete_run.compare --help`` for more details on the
available flags and their usage.

The default baseline is the accepted ``latest-main`` pointer. Promote a reviewed
main result with ``python -m tests.complete_run.baseline promote --run-dir
<run-dir> --channel main``.
"""

from __future__ import annotations

import argparse
import difflib
import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from e3sm_diags.logger import _setup_child_logger, _setup_root_logger
from tests.complete_run.baseline import _load_manifest, _ManifestError
from tests.complete_run.helpers import (
    ComparisonIssue,
    ComparisonSummary,
    compare_netcdf_trees,
    compare_png_trees,
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
DEFAULT_IMAGE_MISMATCH_THRESHOLD = 0.0002

# ``latest-main`` is an accepted-baseline pointer maintained by the promotion
# workflow. It may be a symlink to an immutable complete-run result directory.
DEFAULT_BASELINE_DIR = Path(DEFAULT_RESULTS_DIR) / "latest-main"


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for manual netCDF comparison."""
    _setup_root_logger()

    parser = _build_parser()
    args = parser.parse_args(argv)
    compare_files, compare_values, compare_images = _normalize_modes(args.mode)

    if not compare_files and not compare_values and not compare_images:
        raise ValueError("At least one compare mode must be selected.")

    _validate_compare_dirs(args.dev_dir, args.baseline_dir)
    environment_comparison = _warn_environment_differences(
        args.dev_dir, args.baseline_dir
    )
    report_path = _comparison_report_path(
        args.dev_dir, args.baseline_dir, args.report_dir
    )

    diff_artifact_dir = None
    if args.write_diff_pngs:
        diff_artifact_dir = args.diff_artifact_dir or str(
            report_path.parent / "diff-pngs"
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
    if compare_images:
        image_summary = compare_png_trees(
            dev_root=args.dev_dir,
            baseline_root=args.baseline_dir,
            mismatch_threshold=args.image_mismatch_threshold,
            diff_artifact_dir=diff_artifact_dir,
        )
        _add_image_summary(summary, image_summary)

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
    if "all" in show or "images" in show:
        _render_issue_details("Missing dev images", summary.missing_dev_images)
        _render_issue_details(
            "Missing baseline images", summary.missing_baseline_images
        )
        _render_issue_details("Image mismatches", summary.image_mismatches)

    exit_code = 0 if summary.failure_count == 0 else 1
    _write_comparison_report(
        report_path=report_path,
        dev_dir=args.dev_dir,
        baseline_dir=args.baseline_dir,
        atol=args.atol,
        rtol=args.rtol,
        image_mismatch_threshold=args.image_mismatch_threshold,
        modes=args.mode,
        diff_artifact_dir=diff_artifact_dir,
        environment_comparison=environment_comparison,
        summary=summary,
        exit_code=exit_code,
    )
    logger.info("Wrote comparison report: %s", report_path)

    return exit_code


def _build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for netCDF comparison."""
    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Compare manual complete-run netCDF and PNG outputs against a baseline tree."
        ),
        epilog=(
            "Example: python -m tests.complete_run.compare "
            "--dev-dir </path/to/dev> "
            f"--baseline-dir {DEFAULT_BASELINE_DIR} "
            "--show missing-files --show tolerance-failures"
        ),
    )
    parser.add_argument(
        "--dev-dir",
        required=True,
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
        choices=["all", "files", "data", "images"],
        help=(
            "Comparison mode. Use files for netCDF tree matching, data for shared "
            "netCDF values, images for PNG tree and pixel comparison, or all for "
            "all checks (default: all)."
        ),
    )
    parser.add_argument(
        "--image-mismatch-threshold",
        type=float,
        default=DEFAULT_IMAGE_MISMATCH_THRESHOLD,
        help=(
            "Allowed mismatched-pixel fraction for complete-run PNG comparisons "
            f"(default: {DEFAULT_IMAGE_MISMATCH_THRESHOLD})."
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
            "images",
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
        help=(
            "Directory for optional diff PNGs. Defaults to the comparison "
            "report directory's diff-pngs subdirectory (default: None)."
        ),
    )
    parser.add_argument(
        "--report-dir",
        default=None,
        help=(
            "Directory that receives the JSON report. Defaults to a comparison "
            "directory beside the dev result directory (default: None)."
        ),
    )

    return parser


def _normalize_modes(modes: list[str] | None) -> tuple[bool, bool, bool]:
    """Normalize CLI compare modes into netCDF and image comparison toggles."""
    if not modes or "all" in modes:
        return True, True, True

    compare_files = "files" in modes
    compare_values = "data" in modes

    return compare_files, compare_values, "images" in modes


def _add_image_summary(
    summary: ComparisonSummary, image_summary: ComparisonSummary
) -> None:
    """Add PNG comparison results to an existing netCDF comparison summary."""
    summary.matching_images.extend(image_summary.matching_images)
    summary.missing_dev_images.extend(image_summary.missing_dev_images)
    summary.missing_baseline_images.extend(image_summary.missing_baseline_images)
    summary.image_mismatches.extend(image_summary.image_mismatches)


def _validate_compare_dirs(dev_dir: str | Path, baseline_dir: str | Path) -> None:
    """Validate the directory inputs for the compare workflow.

    Raises
    ------
    FileNotFoundError
        If either comparison root directory does not exist.
    """
    dev_path = Path(dev_dir)
    baseline_path = Path(baseline_dir)

    if baseline_path == DEFAULT_BASELINE_DIR and not baseline_path.is_dir():
        raise FileNotFoundError(
            "No accepted main baseline is promoted: the latest-main baseline "
            f"pointer is missing or broken at {baseline_path}. Review a main "
            "complete-run result, then promote it with: python -m "
            "tests.complete_run.baseline promote --run-dir <run-dir> --channel main"
        )

    missing_dirs = [
        str(path) for path in (dev_path, baseline_path) if not path.is_dir()
    ]

    if missing_dirs:
        missing_message = "\n".join(missing_dirs)
        raise FileNotFoundError(
            f"One or more compare directories do not exist:\n{missing_message}"
        )


def _warn_environment_differences(
    dev_dir: str | Path, baseline_dir: str | Path
) -> dict[str, object]:
    """Log non-blocking dependency provenance differences when available."""
    try:
        dev_manifest = _load_manifest(dev_dir)
        baseline_manifest = _load_manifest(baseline_dir)
    except (FileNotFoundError, _ManifestError) as error:
        logger.info(
            "Environment provenance unavailable; continuing comparison: %s", error
        )
        return {"available": False, "detail": str(error)}

    dev_environment = dev_manifest["environment"]
    baseline_environment = baseline_manifest["environment"]
    differences: list[str] = []
    for key in (
        "python_version",
        "python_implementation",
        "platform",
        "conda_environment",
    ):
        if dev_environment[key] != baseline_environment[key]:
            differences.append(
                f"{key} ({baseline_environment[key]} -> {dev_environment[key]})"
            )
    for package, baseline_version in baseline_environment["packages"].items():
        dev_version = dev_environment["packages"][package]
        if dev_version != baseline_version:
            differences.append(f"{package} ({baseline_version} -> {dev_version})")

    baseline_environment_file = _environment_provenance_path(baseline_dir)
    dev_environment_file = _environment_provenance_path(dev_dir)
    environment_file_diff = _diff_environment_files(
        baseline_environment_file, dev_environment_file
    )
    if differences:
        logger.warning(
            "Complete-run environment provenance differs; review before interpreting "
            "numerical differences: %s\n"
            "baseline environment.yml: %s\n"
            "dev environment.yml: %s\n"
            "environment.yml diff (top-level name excluded):\n%s",
            "; ".join(differences),
            baseline_environment_file,
            dev_environment_file,
            _format_environment_file_diff(
                environment_file_diff,
                baseline_environment_file,
                dev_environment_file,
            ),
        )

    return {
        "available": True,
        "differences": differences,
        "baseline_environment_file": str(baseline_environment_file),
        "dev_environment_file": str(dev_environment_file),
        "environment_file_diff": environment_file_diff,
    }


def _environment_provenance_path(run_dir: str | Path) -> Path:
    """Return the absolute environment.yml provenance path for a complete run."""
    return Path(run_dir).resolve() / "prov" / "environment.yml"


def _diff_environment_files(baseline_path: Path, dev_path: Path) -> dict[str, object]:
    """Return structured environment.yml changes with the top-level name omitted."""
    try:
        baseline_lines = baseline_path.read_text(encoding="utf-8").splitlines()
        dev_lines = dev_path.read_text(encoding="utf-8").splitlines()
    except OSError as error:
        return {"available": False, "detail": str(error), "changes": []}

    baseline_without_name = [
        line for line in baseline_lines if not line.startswith("name:")
    ]
    dev_without_name = [line for line in dev_lines if not line.startswith("name:")]
    matcher = difflib.SequenceMatcher(
        a=baseline_without_name, b=dev_without_name, autojunk=False
    )
    changes = [
        {
            "operation": operation,
            "baseline_lines": baseline_without_name[baseline_start:baseline_end],
            "dev_lines": dev_without_name[dev_start:dev_end],
        }
        for (
            operation,
            baseline_start,
            baseline_end,
            dev_start,
            dev_end,
        ) in matcher.get_opcodes()
        if operation != "equal"
    ]
    return {"available": True, "changes": changes}


def _format_environment_file_diff(
    environment_file_diff: dict[str, object],
    baseline_path: Path,
    dev_path: Path,
) -> str:
    """Render structured environment-file changes as a unified diff for logs."""
    if not environment_file_diff["available"]:
        return f"unavailable: {environment_file_diff['detail']}"

    changes = environment_file_diff["changes"]
    if not changes:
        return "no content differences"

    baseline_lines = baseline_path.read_text(encoding="utf-8").splitlines()
    dev_lines = dev_path.read_text(encoding="utf-8").splitlines()
    baseline_without_name = [
        line for line in baseline_lines if not line.startswith("name:")
    ]
    dev_without_name = [line for line in dev_lines if not line.startswith("name:")]
    diff = difflib.unified_diff(
        baseline_without_name,
        dev_without_name,
        fromfile=str(baseline_path),
        tofile=str(dev_path),
        lineterm="",
    )
    return "\n".join(diff)


def _comparison_report_path(
    dev_dir: str | Path, baseline_dir: str | Path, report_dir: str | Path | None
) -> Path:
    """Build the output path for a comparison's JSON report."""
    root = (
        Path(report_dir).resolve()
        if report_dir is not None
        else Path(dev_dir).resolve().parent / "comparison"
    )
    comparison_name = (
        f"{Path(dev_dir).resolve().name}-vs-{Path(baseline_dir).resolve().name}"
    )
    return root / comparison_name / "comparison-report.json"


def _write_comparison_report(
    *,
    report_path: Path,
    dev_dir: str | Path,
    baseline_dir: str | Path,
    atol: float,
    rtol: float,
    image_mismatch_threshold: float,
    modes: list[str] | None,
    diff_artifact_dir: str | None,
    environment_comparison: dict[str, object],
    summary: ComparisonSummary,
    exit_code: int,
) -> None:
    """Write a JSON record of a complete-run comparison."""
    report = {
        "schema_version": 1,
        "created_at_utc": datetime.now(timezone.utc).isoformat(),
        "status": "passed" if exit_code == 0 else "failed",
        "exit_code": exit_code,
        "paths": {
            "dev_dir": str(Path(dev_dir).resolve()),
            "baseline_dir": str(Path(baseline_dir).resolve()),
            "report": str(report_path),
            "diff_artifact_dir": diff_artifact_dir,
        },
        "comparison_settings": {
            "atol": atol,
            "rtol": rtol,
            "image_mismatch_threshold": image_mismatch_threshold,
            "modes": modes,
        },
        "environment": environment_comparison,
        "summary": {
            "matching_files": [str(path) for path in summary.matching_files],
            "missing_dev_files": [str(path) for path in summary.missing_dev_files],
            "missing_baseline_files": [
                str(path) for path in summary.missing_baseline_files
            ],
            "missing_variables": _issues_to_report(summary.missing_variables),
            "nan_location_mismatches": _issues_to_report(
                summary.nan_location_mismatches
            ),
            "shape_mismatches": _issues_to_report(summary.shape_mismatches),
            "tolerance_failures": _issues_to_report(summary.tolerance_failures),
            "matching_images": [str(path) for path in summary.matching_images],
            "missing_dev_images": [str(path) for path in summary.missing_dev_images],
            "missing_baseline_images": [
                str(path) for path in summary.missing_baseline_images
            ],
            "image_mismatches": _issues_to_report(summary.image_mismatches),
            "compared_file_count": summary.compared_file_count,
            "failure_count": summary.failure_count,
        },
    }
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")


def _issues_to_report(issues: Sequence[ComparisonIssue]) -> list[dict[str, str | None]]:
    """Convert structured comparison issues into JSON-serializable dictionaries."""
    return [
        {
            "relative_path": str(issue.relative_path),
            "var_key": issue.var_key,
            "detail": issue.detail,
            "artifact_path": str(issue.artifact_path)
            if issue.artifact_path is not None
            else None,
        }
        for issue in issues
    ]


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
    logger.info("matched images: %s", len(summary.matching_images))
    logger.info("missing dev images: %s", len(summary.missing_dev_images))
    logger.info("missing baseline images: %s", len(summary.missing_baseline_images))
    logger.info("image mismatches: %s", len(summary.image_mismatches))
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
