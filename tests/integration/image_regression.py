from __future__ import annotations

import json
import os
import subprocess
import sys
from importlib import import_module, metadata
from pathlib import Path
from typing import Any

from tests.integration.utils import _compare_images

BASELINE_METADATA_FILENAME = "baseline_metadata.json"
RUNTIME_METADATA_FILENAME = "runtime_metadata.json"
DEPENDENCY_DIFF_FILENAME = "dependency_diff.json"

_DIST_NAMES_BY_KEY: dict[str, tuple[str, ...]] = {
    "numpy": ("numpy",),
    "pandas": ("pandas",),
    "matplotlib": ("matplotlib",),
    "cartopy": ("cartopy",),
    "xarray": ("xarray",),
    "xcdat": ("xcdat",),
    "xesmf": ("xesmf",),
    "esmf": ("esmf",),
    "esmpy": ("esmpy",),
    "xgcm": ("xgcm",),
}
_MODULE_NAMES_BY_KEY: dict[str, tuple[str, ...]] = {
    "esmf": ("ESMF", "esmpy"),
    "esmpy": ("esmpy", "ESMF"),
}


def assert_image_matches_baseline(
    image_name: str,
    path_to_actual_png: str | Path,
    path_to_expected_png: str | Path,
    baseline_metadata_path: str | Path,
    artifact_dir: str | Path,
    mismatch_threshold: float = 0.0002,
) -> Path:
    """Assert that an image matches its committed baseline.

    Parameters
    ----------
    image_name : str
        Identifier to use for reporting and copied diff artifacts.
    path_to_actual_png : str | Path
        Path to the generated image.
    path_to_expected_png : str | Path
        Path to the expected baseline image.
    baseline_metadata_path : str | Path
        Path to the committed baseline metadata JSON file.
    artifact_dir : str | Path
        Directory used to store diff artifacts and runtime metadata.
    mismatch_threshold : float, optional
        Allowed mismatched pixel fraction before failing, by default 0.0002.

    Returns
    -------
    Path
        Path to the written runtime metadata JSON file.
    """
    actual_path = Path(path_to_actual_png)
    expected_path = Path(path_to_expected_png)
    baseline_metadata = Path(baseline_metadata_path)
    diff_dir = Path(artifact_dir)

    assert actual_path.exists(), f"Generated image does not exist: {actual_path}"
    assert expected_path.exists(), f"Baseline image does not exist: {expected_path}"
    assert baseline_metadata.exists(), (
        f"Baseline metadata does not exist: {baseline_metadata}"
    )

    diff_dir.mkdir(parents=True, exist_ok=True)
    runtime_metadata_path = diff_dir / RUNTIME_METADATA_FILENAME
    write_runtime_metadata(runtime_metadata_path)

    mismatched_images: list[str] = []
    mismatch_fractions: dict[str, float] = {}
    _compare_images(
        mismatched_images,
        image_name,
        str(actual_path),
        str(expected_path),
        diff_dir=str(diff_dir),
        mismatch_threshold=mismatch_threshold,
        mismatch_fractions=mismatch_fractions,
    )
    mismatch_fraction = mismatch_fractions.get(image_name)
    mismatch_fraction_text = (
        f"{mismatch_fraction:.6g}" if mismatch_fraction is not None else "unknown"
    )

    dependency_diff_path = diff_dir / DEPENDENCY_DIFF_FILENAME
    dependency_diff_summary = ""
    if mismatched_images:
        dependency_differences = compare_metadata(
            baseline_metadata,
            runtime_metadata_path,
        )
        write_dependency_diff(
            dependency_diff_path,
            baseline_metadata,
            runtime_metadata_path,
            dependency_differences,
        )
        dependency_diff_summary = format_dependency_diff_summary(
            dependency_differences,
            dependency_diff_path,
        )

    assert not mismatched_images, (
        f"Image regression mismatch for {mismatched_images} "
        f"(fraction={mismatch_fraction_text}, "
        f"threshold={mismatch_threshold}). "
        f"See diff artifacts in {diff_dir}."
        f" {dependency_diff_summary}"
    )

    return runtime_metadata_path


def collect_runtime_metadata() -> dict[str, Any]:
    """Collect resolved runtime metadata for image regression provenance."""
    metadata_dict: dict[str, Any] = {
        "python": sys.version.split()[0],
        "e3sm_diags_git_sha": _get_git_sha(),
        "e3sm_unified_version": os.environ.get("E3SM_UNIFIED_VERSION"),
        "e3sm_unified_feedstock_ref": os.environ.get("E3SM_UNIFIED_FEEDSTOCK_REF"),
        "e3sm_unified_recipe_version": os.environ.get("E3SM_UNIFIED_RECIPE_VERSION"),
        "e3sm_unified_recipe_url": os.environ.get("E3SM_UNIFIED_RECIPE_URL"),
    }

    for key, dist_names in _DIST_NAMES_BY_KEY.items():
        metadata_dict[key] = _get_installed_version(key, dist_names)

    return metadata_dict


def write_runtime_metadata(output_path: str | Path) -> Path:
    """Write runtime dependency metadata to a JSON file."""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(collect_runtime_metadata(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return path


def compare_metadata(
    baseline_metadata_path: str | Path,
    runtime_metadata_path: str | Path,
) -> dict[str, dict[str, Any]]:
    baseline_metadata = json.loads(
        Path(baseline_metadata_path).read_text(encoding="utf-8")
    )
    runtime_metadata = json.loads(
        Path(runtime_metadata_path).read_text(encoding="utf-8")
    )

    differences: dict[str, dict[str, Any]] = {}
    for key in sorted(set(baseline_metadata) | set(runtime_metadata)):
        baseline_value = baseline_metadata.get(key)
        runtime_value = runtime_metadata.get(key)
        if baseline_value != runtime_value:
            differences[key] = {
                "baseline": baseline_value,
                "runtime": runtime_value,
            }

    return differences


def write_dependency_diff(
    output_path: str | Path,
    baseline_metadata_path: str | Path,
    runtime_metadata_path: str | Path,
    differences: dict[str, dict[str, Any]],
) -> Path:
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "baseline_metadata_path": str(baseline_metadata_path),
        "runtime_metadata_path": str(runtime_metadata_path),
        "differences": differences,
    }
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    return path


def format_dependency_diff_summary(
    differences: dict[str, dict[str, Any]],
    dependency_diff_path: str | Path,
) -> str:
    if not differences:
        return (
            "No dependency metadata differences were detected. "
            f"See {dependency_diff_path}."
        )

    changed_dependencies = "; ".join(
        f"{key}: {values['baseline']} -> {values['runtime']}"
        for key, values in differences.items()
    )
    return (
        f"Dependency differences: {changed_dependencies}. See {dependency_diff_path}."
    )


def _get_git_sha() -> str | None:
    try:
        proc = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None

    return proc.stdout.strip() or None


def _get_installed_version(key: str, dist_names: tuple[str, ...]) -> str | None:
    for dist_name in dist_names:
        try:
            return metadata.version(dist_name)
        except metadata.PackageNotFoundError:
            continue

    for module_name in _MODULE_NAMES_BY_KEY.get(key, ()):
        try:
            module = import_module(module_name)
        except ModuleNotFoundError:
            continue

        version = getattr(module, "__version__", getattr(module, "VERSION", None))
        if version is None:
            continue

        if isinstance(version, tuple):
            return ".".join(str(item) for item in version)

        return str(version)

    return None


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("output", type=Path)
    args = parser.parse_args()

    write_runtime_metadata(args.output)
