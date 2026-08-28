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
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, Literal

import numpy as np
from PIL import Image, ImageChops, ImageDraw

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
    matching_images: list[Path] = field(default_factory=list)
    missing_dev_images: list[Path] = field(default_factory=list)
    missing_baseline_images: list[Path] = field(default_factory=list)
    image_mismatches: list[ComparisonIssue] = field(default_factory=list)

    @property
    def compared_file_count(self) -> int:
        """Total number of shared files that were value-compared."""
        failed_paths = {
            issue.relative_path
            for issues in (
                self.missing_variables,
                self.nan_location_mismatches,
                self.shape_mismatches,
                self.tolerance_failures,
            )
            for issue in issues
        }
        return len(self.matching_files) + len(failed_paths)

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
            + len(self.missing_dev_images)
            + len(self.missing_baseline_images)
            + len(self.image_mismatches)
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


def _discover_png_files(root_dir: str | Path) -> dict[Path, Path]:
    """Discover diagnostic PNG files below a complete-run root directory."""
    root_path = Path(root_dir).expanduser().resolve()

    return {
        file_path.relative_to(root_path): file_path
        for file_path in sorted(root_path.rglob("*.png"))
    }


def match_png_files(dev_root: str | Path, baseline_root: str | Path) -> FileTreeMatch:
    """Match diagnostic PNG files by relative path across two directory trees."""
    dev_paths = set(_discover_png_files(dev_root))
    baseline_paths = set(_discover_png_files(baseline_root))

    return FileTreeMatch(
        shared_paths=sorted(dev_paths & baseline_paths),
        missing_dev_paths=sorted(baseline_paths - dev_paths),
        missing_baseline_paths=sorted(dev_paths - baseline_paths),
    )


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
) -> list[ComparisonOutcome]:
    """Compare every data variable in a shared netCDF file pair."""
    import xarray as xr

    relative = Path(relative_path)

    logger.info("Comparing %s", relative)

    with (
        xr.open_dataset(dev_path) as dev_ds,
        xr.open_dataset(baseline_path) as baseline_ds,
    ):
        return compare_dataset_pair(
            dev_ds,
            baseline_ds,
            relative_path=relative,
            atol=atol,
            rtol=rtol,
            diff_artifact_dir=diff_artifact_dir,
        )


def compare_dataset_pair(
    dev_ds: xr.Dataset,
    baseline_ds: xr.Dataset,
    *,
    relative_path: str | Path,
    atol: float,
    rtol: float,
    diff_artifact_dir: str | Path | None = None,
) -> list[ComparisonOutcome]:
    """Compare all data variables shared by two opened datasets."""
    relative = Path(relative_path)
    dev_var_keys = set(dev_ds.data_vars)
    baseline_var_keys = set(baseline_ds.data_vars)
    outcomes: list[ComparisonOutcome] = []

    for var_key in sorted(dev_var_keys - baseline_var_keys):
        outcomes.append(
            ComparisonOutcome(
                status="missing_variable",
                relative_path=relative,
                var_key=var_key,
                detail=f"Variable {var_key!r} is missing from the baseline dataset.",
            )
        )

    for var_key in sorted(baseline_var_keys - dev_var_keys):
        outcomes.append(
            ComparisonOutcome(
                status="missing_variable",
                relative_path=relative,
                var_key=var_key,
                detail=f"Variable {var_key!r} is missing from the dev dataset.",
            )
        )

    for var_key in sorted(dev_var_keys & baseline_var_keys):
        dev_data = dev_ds[var_key].values
        baseline_data = baseline_ds[var_key].values
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
                title=f"{relative} [{var_key}] ({status})",
            )
        outcomes.append(
            ComparisonOutcome(
                status=status,
                relative_path=relative,
                var_key=var_key,
                detail=comparison_detail,
                artifact_path=artifact_path,
            )
        )

    if not outcomes:
        outcomes.append(
            ComparisonOutcome(
                status="missing_variable",
                relative_path=relative,
                detail="Neither dataset contains data variables.",
            )
        )

    return outcomes


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

    # QBO output names predate the standard plot filename convention.  The
    # diagnostic always writes the zonal-wind variable under one of these
    # names, so its variable key cannot be obtained from hyphen-separated
    # filename components.
    if path.stem in {"qbo_diags_qbo_test", "qbo_diags_qbo_ref"}:
        return "U"

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
        outcomes = compare_file_pair(
            dev_files[relative_path],
            baseline_files[relative_path],
            relative_path=relative_path,
            atol=atol,
            rtol=rtol,
            diff_artifact_dir=diff_artifact_dir,
        )

        if all(outcome.status == "matching" for outcome in outcomes):
            summary.matching_files.append(relative_path)
            continue

        for outcome in outcomes:
            if outcome.status == "matching":
                continue
            issue = ComparisonIssue(
                relative_path=outcome.relative_path,
                var_key=outcome.var_key,
                detail=outcome.detail,
                artifact_path=outcome.artifact_path,
            )
            if outcome.status == "missing_variable":
                summary.missing_variables.append(issue)
            elif outcome.status == "nan_location_mismatch":
                summary.nan_location_mismatches.append(issue)
            elif outcome.status == "shape_mismatch":
                summary.shape_mismatches.append(issue)
            else:
                summary.tolerance_failures.append(issue)

    return summary


def compare_png_trees(
    dev_root: str | Path,
    baseline_root: str | Path,
    *,
    mismatch_threshold: float,
    diff_artifact_dir: str | Path | None = None,
) -> ComparisonSummary:
    """Compare complete-run PNG trees using pixel-level image regression.

    This intentionally uses the same mismatched-pixel fraction and default
    threshold as the committed targeted image-regression suite.  A separate
    summary is returned so callers can combine PNG and netCDF validation.
    """
    dev_images = _discover_png_files(dev_root)
    baseline_images = _discover_png_files(baseline_root)
    tree_match = match_png_files(dev_root, baseline_root)
    summary = ComparisonSummary(
        missing_dev_images=tree_match.missing_dev_paths.copy(),
        missing_baseline_images=tree_match.missing_baseline_paths.copy(),
    )

    for relative_path in tree_match.shared_paths:
        mismatch = compare_png_pair(
            dev_images[relative_path],
            baseline_images[relative_path],
            relative_path=relative_path,
            mismatch_threshold=mismatch_threshold,
            diff_artifact_dir=diff_artifact_dir,
        )
        if mismatch is None:
            summary.matching_images.append(relative_path)
        else:
            summary.image_mismatches.append(mismatch)

    return summary


def compare_png_pair(
    dev_path: str | Path,
    baseline_path: str | Path,
    *,
    relative_path: str | Path,
    mismatch_threshold: float,
    diff_artifact_dir: str | Path | None = None,
) -> ComparisonIssue | None:
    """Return an issue when a PNG pair exceeds the mismatch threshold."""
    relative = Path(relative_path)
    with Image.open(dev_path) as dev_image, Image.open(baseline_path) as baseline_image:
        dev_rgb = dev_image.convert("RGB")
        baseline_rgb = baseline_image.convert("RGB")

    if dev_rgb.size != baseline_rgb.size:
        detail = f"Image size mismatch: {dev_rgb.size} != {baseline_rgb.size}."
        return _image_mismatch_issue(
            dev_path, baseline_path, relative, detail, None, diff_artifact_dir
        )

    diff = ImageChops.difference(dev_rgb, baseline_rgb)
    bbox = diff.getbbox()
    if bbox is None:
        return None

    nonzero_pixels = (
        diff.crop(bbox)
        .point(lambda value: 255 if value else 0)
        .convert("L")
        .point(bool)
        .getdata()
    )
    mismatch_fraction = sum(nonzero_pixels) / (baseline_rgb.width * baseline_rgb.height)
    if mismatch_fraction < mismatch_threshold:
        return None

    detail = (
        f"Mismatched pixel fraction: {mismatch_fraction:.6g} "
        f"(threshold: {mismatch_threshold})."
    )
    return _image_mismatch_issue(
        dev_path,
        baseline_path,
        relative,
        detail,
        diff,
        diff_artifact_dir,
    )


def _image_mismatch_issue(
    dev_path: str | Path,
    baseline_path: str | Path,
    relative_path: Path,
    detail: str,
    diff: Image.Image | None,
    diff_artifact_dir: str | Path | None,
) -> ComparisonIssue:
    """Create a PNG mismatch issue and, when requested, its review artifacts."""
    artifact_path = None
    if diff_artifact_dir is not None:
        output_dir = Path(diff_artifact_dir) / "image-diffs" / relative_path.parent
        output_dir.mkdir(parents=True, exist_ok=True)
        stem = relative_path.stem
        shutil.copy(dev_path, output_dir / f"{stem}_actual.png")
        shutil.copy(baseline_path, output_dir / f"{stem}_expected.png")
        artifact_path = output_dir / f"{stem}_diff.png"
        if diff is not None:
            draw = ImageDraw.Draw(diff)
            bbox = diff.getbbox()
            if bbox is not None:
                draw.rectangle(bbox, outline="red")
            diff.save(artifact_path, "PNG")

    return ComparisonIssue(
        relative_path=relative_path,
        detail=detail,
        artifact_path=artifact_path,
    )


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
