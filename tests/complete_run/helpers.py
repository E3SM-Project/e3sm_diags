"""Provide reusable helpers for manual complete-run comparison workflows.

This module contains the pure and reusable building blocks behind the manual
netCDF comparison flow, including file discovery, filename-to-variable
resolution, array-difference classification, comparison summaries, and
optional diff artifact generation. It exists so comparison behavior can be
shared by the CLI module and covered by focused unit tests without requiring
an end-to-end HPC run.
"""

from __future__ import annotations

import re
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
from e3sm_diags.logger import _setup_child_logger

if TYPE_CHECKING:
    import matplotlib.figure
    import xarray as xr

logger = _setup_child_logger(__name__)

REPO_ROOT = Path(__file__).resolve().parents[2]

ComparisonStatus = Literal[
    "matching",
    "missing_variable",
    "nan_location_mismatch",
    "shape_mismatch",
    "tolerance_failure",
]


@dataclass(frozen=True)
class ComparisonIssue:
    """A single comparison outcome for a netCDF file."""

    relative_path: Path
    var_key: str | None = None
    detail: str | None = None
    artifact_path: Path | None = None


@dataclass(frozen=True)
class ComparisonOutcome:
    """Structured result for comparing one netCDF file pair."""

    status: ComparisonStatus
    relative_path: Path
    var_key: str | None = None
    detail: str | None = None
    artifact_path: Path | None = None


@dataclass(frozen=True)
class FileTreeMatch:
    """Relative-path alignment between two netCDF directory trees."""

    shared_paths: list[Path]
    missing_dev_paths: list[Path]
    missing_baseline_paths: list[Path]


@dataclass
class ComparisonSummary:
    """Aggregate results across a directory-tree comparison."""

    matching_files: list[Path] = field(default_factory=list)
    missing_dev_files: list[Path] = field(default_factory=list)
    missing_baseline_files: list[Path] = field(default_factory=list)
    missing_variables: list[ComparisonIssue] = field(default_factory=list)
    nan_location_mismatches: list[ComparisonIssue] = field(default_factory=list)
    shape_mismatches: list[ComparisonIssue] = field(default_factory=list)
    tolerance_failures: list[ComparisonIssue] = field(default_factory=list)

    @property
    def compared_file_count(self) -> int:
        """Total number of shared files that were value-compared."""
        return (
            len(self.matching_files)
            + len(self.missing_variables)
            + len(self.nan_location_mismatches)
            + len(self.shape_mismatches)
            + len(self.tolerance_failures)
        )

    @property
    def failure_count(self) -> int:
        """Total number of comparison failures."""
        return (
            len(self.missing_dev_files)
            + len(self.missing_baseline_files)
            + len(self.missing_variables)
            + len(self.nan_location_mismatches)
            + len(self.shape_mismatches)
            + len(self.tolerance_failures)
        )


def append_run_suffix(
    results_dir: str | Path,
) -> str:
    """Append a generated run suffix to a complete-run results root."""
    return f"{Path(results_dir).as_posix().rstrip('/')}/{_build_run_suffix()}"


def _build_run_suffix() -> str:
    """Build the timestamp-branch-hash suffix for complete-run outputs."""
    timestamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    branch_name = _sanitize_results_dir_component(
        _get_git_output(["rev-parse", "--abbrev-ref", "HEAD"], "unknown-branch")
    )
    commit_hash = _sanitize_results_dir_component(
        _get_git_output(["rev-parse", "--short", "HEAD"], "unknown-hash")
    )

    return f"{timestamp}-{branch_name}-{commit_hash}/"


def _sanitize_results_dir_component(value: str) -> str:
    sanitized = re.sub(r"[^A-Za-z0-9._-]+", "-", value).strip("-")

    return sanitized or "unknown"


def _get_git_output(args: list[str], fallback: str) -> str:
    try:
        completed = subprocess.run(
            ["git", *args],
            check=True,
            capture_output=True,
            text=True,
            cwd=REPO_ROOT,
        )
    except (subprocess.CalledProcessError, FileNotFoundError):
        return fallback

    output = completed.stdout.strip()

    return output or fallback


def match_netcdf_files(
    dev_root: str | Path,
    baseline_root: str | Path,
) -> FileTreeMatch:
    """Match netCDF files by relative path across two directory trees."""
    dev_files = _discover_netcdf_files(dev_root)
    baseline_files = _discover_netcdf_files(baseline_root)

    dev_paths = set(dev_files)
    baseline_paths = set(baseline_files)

    return FileTreeMatch(
        shared_paths=sorted(dev_paths & baseline_paths),
        missing_dev_paths=sorted(baseline_paths - dev_paths),
        missing_baseline_paths=sorted(dev_paths - baseline_paths),
    )


def _discover_netcdf_files(root_dir: str | Path) -> dict[Path, Path]:
    """Discover netCDF files below a root directory.

    Parameters
    ----------
    root_dir : str | Path
        Root directory containing comparison outputs.

    Returns
    -------
    dict[Path, Path]
        Mapping of relative paths to absolute file paths.
    """
    root_path = Path(root_dir).expanduser().resolve()

    return {
        file_path.relative_to(root_path): file_path
        for file_path in sorted(root_path.rglob("*.nc"))
    }


def expand_candidate_var_keys(var_key: str) -> list[str]:
    """Expand a target variable into raw and derived candidate names."""
    candidates: list[str] = [var_key]

    derived_map = DERIVED_VARIABLES.get(var_key) or DERIVED_VARIABLES.get(
        var_key.upper()
    )
    if derived_map is not None:
        for source_keys in derived_map:
            for source_key in source_keys:
                if source_key not in candidates:
                    candidates.append(source_key)

    return candidates


def get_var_data(
    ds: xr.Dataset,
    var_key: str,
) -> tuple[np.ndarray | None, str | None]:
    """Retrieve variable data using the current complete-run lookup strategy.

    Parameters
    ----------
    ds : xr.Dataset
        Dataset to inspect.
    var_key : str
        Variable key inferred from the filename.

    Returns
    -------
    tuple[np.ndarray | None, str | None]
        The matched data array values and the variable name used.
    """
    for candidate_key in expand_candidate_var_keys(var_key):
        if candidate_key in ds.data_vars:
            return ds[candidate_key].values, candidate_key

    return None, None


def classify_array_difference(
    dev_data: np.ndarray,
    baseline_data: np.ndarray,
    *,
    atol: float,
    rtol: float,
) -> tuple[ComparisonStatus, str | None]:
    """Classify two arrays according to the current manual compare workflow."""
    if dev_data.shape != baseline_data.shape:
        return (
            "shape_mismatch",
            f"Shape mismatch: {dev_data.shape} != {baseline_data.shape}",
        )

    dev_nan_mask = np.isnan(dev_data)
    baseline_nan_mask = np.isnan(baseline_data)
    if not np.array_equal(dev_nan_mask, baseline_nan_mask):
        mismatch_count = int(np.count_nonzero(dev_nan_mask != baseline_nan_mask))
        return (
            "nan_location_mismatch",
            f"NaN locations differ at {mismatch_count} element(s).",
        )

    try:
        np.testing.assert_allclose(
            dev_data,
            baseline_data,
            atol=atol,
            rtol=rtol,
            equal_nan=True,
        )
    except AssertionError as err:
        return ("tolerance_failure", str(err))

    return ("matching", None)


def compare_file_pair(
    dev_path: str | Path,
    baseline_path: str | Path,
    *,
    relative_path: str | Path,
    atol: float,
    rtol: float,
    diff_artifact_dir: str | Path | None = None,
) -> ComparisonOutcome:
    """Compare a shared netCDF file pair."""
    import xarray as xr

    relative = Path(relative_path)
    var_key = infer_variable_key_from_path(relative)

    logger.info("Comparing %s", relative)

    with (
        xr.open_dataset(dev_path) as dev_ds,
        xr.open_dataset(baseline_path) as baseline_ds,
    ):
        dev_data, dev_match = get_var_data(dev_ds, var_key)
        baseline_data, baseline_match = get_var_data(baseline_ds, var_key)

        if dev_data is None or baseline_data is None:
            missing_side = "dev" if dev_data is None else "baseline"
            missing_detail = (
                f"Missing variable for {missing_side} dataset. "
                f"Requested {var_key!r}, matched dev={dev_match!r}, baseline={baseline_match!r}."
            )
            return ComparisonOutcome(
                status="missing_variable",
                relative_path=relative,
                var_key=var_key,
                detail=missing_detail,
            )

        status, comparison_detail = classify_array_difference(
            dev_data,
            baseline_data,
            atol=atol,
            rtol=rtol,
        )
        artifact_path = None
        if status != "matching" and diff_artifact_dir is not None:
            artifact_path = write_diff_artifact(
                dev_data=dev_data,
                baseline_data=baseline_data,
                artifact_root=diff_artifact_dir,
                relative_path=relative,
                var_key=var_key,
                title=f"{relative} ({status})",
            )

        return ComparisonOutcome(
            status=status,
            relative_path=relative,
            var_key=var_key,
            detail=comparison_detail,
            artifact_path=artifact_path,
        )


def infer_variable_key_from_path(relative_path: str | Path) -> str:
    """Infer the primary variable key from a diagnostic netCDF filename.

    Parameters
    ----------
    relative_path : str | Path
        Relative file path for a netCDF artifact.

    Returns
    -------
    str
        The inferred variable key.

    Raises
    ------
    ValueError
        If the filename does not contain enough ``-``-separated parts.
    """
    path = Path(relative_path)
    parts = path.stem.split("-")

    if len(parts) < 3:
        raise ValueError(f"Could not infer variable key from filename: {path.name!r}.")

    var_key = parts[-3]
    if var_key.isdigit():
        if len(parts) < 4:
            raise ValueError(
                f"Could not infer 3D variable key from filename: {path.name!r}."
            )
        var_key = parts[-4]

    return var_key


def compare_netcdf_trees(
    dev_root: str | Path,
    baseline_root: str | Path,
    *,
    atol: float,
    rtol: float,
    compare_values: bool = True,
    diff_artifact_dir: str | Path | None = None,
) -> ComparisonSummary:
    """Compare netCDF trees and collect a summary of the results."""
    dev_files = _discover_netcdf_files(dev_root)
    baseline_files = _discover_netcdf_files(baseline_root)
    tree_match = match_netcdf_files(dev_root, baseline_root)

    summary = ComparisonSummary(
        missing_dev_files=tree_match.missing_dev_paths.copy(),
        missing_baseline_files=tree_match.missing_baseline_paths.copy(),
    )

    if not compare_values:
        return summary

    for relative_path in tree_match.shared_paths:
        outcome = compare_file_pair(
            dev_files[relative_path],
            baseline_files[relative_path],
            relative_path=relative_path,
            atol=atol,
            rtol=rtol,
            diff_artifact_dir=diff_artifact_dir,
        )

        if outcome.status == "matching":
            summary.matching_files.append(outcome.relative_path)
        elif outcome.status == "missing_variable":
            summary.missing_variables.append(
                ComparisonIssue(
                    relative_path=outcome.relative_path,
                    var_key=outcome.var_key,
                    detail=outcome.detail,
                    artifact_path=outcome.artifact_path,
                )
            )
        elif outcome.status == "nan_location_mismatch":
            summary.nan_location_mismatches.append(
                ComparisonIssue(
                    relative_path=outcome.relative_path,
                    var_key=outcome.var_key,
                    detail=outcome.detail,
                    artifact_path=outcome.artifact_path,
                )
            )
        elif outcome.status == "shape_mismatch":
            summary.shape_mismatches.append(
                ComparisonIssue(
                    relative_path=outcome.relative_path,
                    var_key=outcome.var_key,
                    detail=outcome.detail,
                    artifact_path=outcome.artifact_path,
                )
            )
        else:
            summary.tolerance_failures.append(
                ComparisonIssue(
                    relative_path=outcome.relative_path,
                    var_key=outcome.var_key,
                    detail=outcome.detail,
                    artifact_path=outcome.artifact_path,
                )
            )

    return summary


def write_diff_artifact(
    *,
    dev_data: np.ndarray,
    baseline_data: np.ndarray,
    artifact_root: str | Path,
    relative_path: str | Path,
    var_key: str,
    title: str,
) -> Path:
    """Write a PNG artifact that visualizes the dev-baseline differences."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    artifact_root_path = Path(artifact_root)
    relative = Path(relative_path)
    output_dir = artifact_root_path / relative.parent
    output_dir.mkdir(parents=True, exist_ok=True)

    output_path = output_dir / f"{relative.stem}-{var_key}-diff.png"
    figure = _build_diff_figure(dev_data, baseline_data, title=title)
    figure.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(figure)

    return output_path


def _build_diff_figure(
    dev_data: np.ndarray,
    baseline_data: np.ndarray,
    *,
    title: str,
) -> matplotlib.figure.Figure:
    """Build a simple PNG figure for mismatch debugging."""
    import matplotlib.pyplot as plt

    dev_flat = np.ravel(dev_data.astype(float, copy=False))
    baseline_flat = np.ravel(baseline_data.astype(float, copy=False))
    diff_flat = dev_flat - baseline_flat

    figure, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)
    axes[0].plot(dev_flat, linewidth=0.8)
    axes[0].set_ylabel("dev")
    axes[1].plot(baseline_flat, linewidth=0.8)
    axes[1].set_ylabel("baseline")
    axes[2].plot(diff_flat, linewidth=0.8)
    axes[2].set_ylabel("diff")
    axes[2].set_xlabel("flattened index")
    figure.suptitle(title)

    return figure
