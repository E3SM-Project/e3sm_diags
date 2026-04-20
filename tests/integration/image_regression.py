from __future__ import annotations

import json
import os
import subprocess
import sys
from importlib import import_module
from importlib import metadata
from pathlib import Path
from typing import Any

from tests.integration.utils import _compare_images

BASELINE_METADATA_FILENAME = "baseline_metadata.json"
RUNTIME_METADATA_FILENAME = "runtime_metadata.json"

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
    diff_dir = Path(artifact_dir)

    assert actual_path.exists(), f"Generated image does not exist: {actual_path}"
    assert expected_path.exists(), f"Baseline image does not exist: {expected_path}"

    diff_dir.mkdir(parents=True, exist_ok=True)
    runtime_metadata_path = diff_dir / RUNTIME_METADATA_FILENAME
    write_runtime_metadata(runtime_metadata_path)

    mismatched_images: list[str] = []
    _compare_images(
        mismatched_images,
        image_name,
        str(actual_path),
        str(expected_path),
        diff_dir=str(diff_dir),
        mismatch_threshold=mismatch_threshold,
    )

    assert len(mismatched_images) == 0

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
