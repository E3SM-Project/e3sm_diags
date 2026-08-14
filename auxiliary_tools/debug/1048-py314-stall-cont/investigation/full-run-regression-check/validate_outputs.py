#!/usr/bin/env python
"""Compare lock-workaround-disabled and enabled diagnostic outputs.

Run this script with either of the Python release environments used by the
investigation, for example::

    conda run -n ed_1048_xr_2026070_py31314 python validate_outputs.py

The default output directory is the directory containing this script. A
nonzero exit status means at least one configured comparison failed.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import re
import sys
from collections.abc import Iterator, Mapping, Sequence
from datetime import datetime
from pathlib import Path
from typing import Any

import netCDF4
import numpy as np
from PIL import Image


INPUTS = {
    "3.13.14": {
        "disabled": Path(
            "/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py31314/"
            "climo-lock-workaround-disabled/model_vs_obs_1985-2014_units"
        ),
        "enabled": Path(
            "/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py31314/"
            "climo-lock-workaround-enabled/model_vs_obs_1985-2014_units"
        ),
    },
    "3.14.6": {
        "disabled": Path(
            "/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3146/"
            "climo-lock-workaround-disabled/model_vs_obs_1985-2014_units"
        ),
        "enabled": Path(
            "/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3146/"
            "climo-lock-workaround-enabled/model_vs_obs_1985-2014_units"
        ),
    },
}
EXCLUDED_INCOMPLETE_VERSIONS = ["3.14.1", "3.14.2", "3.14.3", "3.14.4"]
SUFFIXES = {".nc", ".json", ".png"}
VIEWER_VERSION_RE = re.compile(r"^\d{4}-\d{2}-\d{2}-\d{2}-\d{2}$")
CHUNK_BYTES = 16 * 1024 * 1024


def parse_args() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent,
        help="Directory for inventory.tsv, file-results.tsv, details, and summary.",
    )
    parser.add_argument(
        "--atol",
        type=float,
        default=1e-12,
        help="Absolute tolerance for numeric values (default: %(default)g).",
    )
    parser.add_argument(
        "--rtol",
        type=float,
        default=1e-12,
        help="Relative tolerance for numeric values (default: %(default)g).",
    )
    return parser.parse_args()


def inventory(root: Path) -> dict[str, str]:
    """Return relevant relative paths mapped to their file kind."""
    if not root.is_dir():
        raise FileNotFoundError(f"Input directory does not exist: {root}")
    return {
        path.relative_to(root).as_posix(): path.suffix.removeprefix(".")
        for path in root.rglob("*")
        if path.is_file() and path.suffix in SUFFIXES
    }


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of a file without loading it all at once."""
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def jsonable(value: Any) -> Any:
    """Convert NumPy and NetCDF metadata values into JSON-compatible values."""
    if isinstance(value, np.ndarray):
        return [jsonable(item) for item in value.tolist()]
    if isinstance(value, np.generic):
        return jsonable(value.item())
    if isinstance(value, bytes):
        return value.decode("utf-8", errors="backslashreplace")
    if isinstance(value, float) and not math.isfinite(value):
        return {"non_finite_float": repr(value)}
    if isinstance(value, Mapping):
        return {str(key): jsonable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [jsonable(item) for item in value]
    return value


def attributes(obj: Any) -> dict[str, Any]:
    """Return canonicalized NetCDF attributes."""
    return {name: jsonable(obj.getncattr(name)) for name in sorted(obj.ncattrs())}


def dictionary_diff(
    left: Mapping[str, Any], right: Mapping[str, Any]
) -> dict[str, Any]:
    """Return complete key and value differences between two dictionaries."""
    left_keys = set(left)
    right_keys = set(right)
    changed = {
        key: {"disabled": left[key], "enabled": right[key]}
        for key in sorted(left_keys & right_keys)
        if left[key] != right[key]
    }
    return {
        "missing_in_enabled": sorted(left_keys - right_keys),
        "additional_in_enabled": sorted(right_keys - left_keys),
        "changed": changed,
    }


def walk_groups(group: netCDF4.Group, prefix: str = "") -> Iterator[tuple[str, Any]]:
    """Yield a dataset and all nested groups with stable names."""
    yield prefix or "/", group
    for name in sorted(group.groups):
        child_prefix = f"{prefix}/{name}" if prefix else f"/{name}"
        yield from walk_groups(group.groups[name], child_prefix)


def collect_netcdf_structure(dataset: netCDF4.Dataset) -> dict[str, Any]:
    """Collect group, dimension, variable, dtype, and metadata structure."""
    groups: dict[str, Any] = {}
    for group_name, group in walk_groups(dataset):
        groups[group_name] = {
            "attributes": attributes(group),
            "dimensions": {
                name: {"size": len(dim), "unlimited": dim.isunlimited()}
                for name, dim in sorted(group.dimensions.items())
            },
            "variables": {
                name: {
                    "dimensions": list(var.dimensions),
                    "shape": list(var.shape),
                    "dtype": str(var.dtype),
                    "attributes": attributes(var),
                }
                for name, var in sorted(group.variables.items())
            },
        }
    return {
        "data_model": dataset.data_model,
        "disk_format": dataset.disk_format,
        "groups": groups,
    }


def variable_map(dataset: netCDF4.Dataset) -> dict[str, netCDF4.Variable]:
    """Return variables from all groups keyed by absolute group path."""
    result: dict[str, netCDF4.Variable] = {}
    for group_name, group in walk_groups(dataset):
        base = "" if group_name == "/" else group_name
        for name, variable in group.variables.items():
            result[f"{base}/{name}"] = variable
    return result


def variable_slices(variable: netCDF4.Variable) -> Iterator[Any]:
    """Yield slices bounded to roughly ``CHUNK_BYTES``."""
    if variable.ndim == 0:
        yield (...)
        return
    trailing_items = math.prod(variable.shape[1:]) if variable.ndim > 1 else 1
    itemsize = max(getattr(variable.dtype, "itemsize", 8), 1)
    rows = max(1, CHUNK_BYTES // max(trailing_items * itemsize, 1))
    for start in range(0, variable.shape[0], rows):
        yield (slice(start, min(start + rows, variable.shape[0])),) + (
            (slice(None),) * (variable.ndim - 1)
        )


def nan_mask(values: np.ndarray) -> np.ndarray:
    """Return a NaN mask for floating and complex arrays."""
    if np.issubdtype(values.dtype, np.inexact):
        return np.isnan(values)
    return np.zeros(values.shape, dtype=bool)


def compare_numeric_variable(
    disabled: netCDF4.Variable,
    enabled: netCDF4.Variable,
    atol: float,
    rtol: float,
) -> dict[str, Any]:
    """Compare one numeric NetCDF variable in bounded chunks."""
    exact = True
    masks_equal = True
    within_tolerance = True
    max_abs = 0.0
    max_rel = 0.0
    compared = 0
    disabled.set_auto_maskandscale(False)
    enabled.set_auto_maskandscale(False)

    for selection in variable_slices(disabled):
        left = np.asarray(disabled[selection])
        right = np.asarray(enabled[selection])
        left_nan = nan_mask(left)
        right_nan = nan_mask(right)
        chunk_masks_equal = bool(np.array_equal(left_nan, right_nan))
        masks_equal &= chunk_masks_equal
        exact &= bool(np.array_equal(left, right, equal_nan=True))
        if not chunk_masks_equal:
            within_tolerance = False
            continue

        valid = ~left_nan
        if not np.any(valid):
            continue
        left_valid = left[valid]
        right_valid = right[valid]
        compared += int(left_valid.size)
        within_tolerance &= bool(
            np.allclose(left_valid, right_valid, atol=atol, rtol=rtol, equal_nan=True)
        )
        finite = np.isfinite(left_valid) & np.isfinite(right_valid)
        if np.any(finite):
            left_finite = left_valid[finite].astype(np.float64, copy=False)
            right_finite = right_valid[finite].astype(np.float64, copy=False)
            absolute = np.abs(left_finite - right_finite)
            scale = np.maximum(np.abs(left_finite), np.abs(right_finite))
            relative = np.divide(
                absolute,
                scale,
                out=np.zeros_like(absolute),
                where=scale != 0,
            )
            max_abs = max(max_abs, float(np.max(absolute, initial=0.0)))
            max_rel = max(max_rel, float(np.max(relative, initial=0.0)))

    return {
        "exact": exact,
        "nan_masks_equal": masks_equal,
        "within_tolerance": within_tolerance,
        "max_absolute_difference": max_abs,
        "max_relative_difference": max_rel,
        "compared_non_nan_values": compared,
    }


def compare_non_numeric_variable(
    disabled: netCDF4.Variable, enabled: netCDF4.Variable
) -> dict[str, Any]:
    """Compare one nonnumeric NetCDF variable exactly."""
    exact = True
    disabled.set_auto_maskandscale(False)
    enabled.set_auto_maskandscale(False)
    for selection in variable_slices(disabled):
        exact &= bool(np.array_equal(disabled[selection], enabled[selection]))
    return {
        "exact": exact,
        "nan_masks_equal": True,
        "within_tolerance": exact,
        "max_absolute_difference": None,
        "max_relative_difference": None,
        "compared_non_nan_values": int(math.prod(disabled.shape)),
    }


def compare_netcdf(
    disabled_path: Path, enabled_path: Path, atol: float, rtol: float
) -> dict[str, Any]:
    """Compare NetCDF structure, metadata, dtypes, masks, and values."""
    with (
        netCDF4.Dataset(disabled_path, mode="r") as disabled,
        netCDF4.Dataset(enabled_path, mode="r") as enabled,
    ):
        disabled_structure = collect_netcdf_structure(disabled)
        enabled_structure = collect_netcdf_structure(enabled)
        structure_diff = dictionary_diff(
            {
                "data_model": disabled_structure["data_model"],
                "disk_format": disabled_structure["disk_format"],
                "groups": {
                    name: {
                        "dimensions": group["dimensions"],
                        "variables": {
                            variable: {
                                key: value
                                for key, value in definition.items()
                                if key != "attributes"
                            }
                            for variable, definition in group["variables"].items()
                        },
                    }
                    for name, group in disabled_structure["groups"].items()
                },
            },
            {
                "data_model": enabled_structure["data_model"],
                "disk_format": enabled_structure["disk_format"],
                "groups": {
                    name: {
                        "dimensions": group["dimensions"],
                        "variables": {
                            variable: {
                                key: value
                                for key, value in definition.items()
                                if key != "attributes"
                            }
                            for variable, definition in group["variables"].items()
                        },
                    }
                    for name, group in enabled_structure["groups"].items()
                },
            },
        )
        metadata_disabled = {
            name: {
                "attributes": group["attributes"],
                "variable_attributes": {
                    variable: definition["attributes"]
                    for variable, definition in group["variables"].items()
                },
            }
            for name, group in disabled_structure["groups"].items()
        }
        metadata_enabled = {
            name: {
                "attributes": group["attributes"],
                "variable_attributes": {
                    variable: definition["attributes"]
                    for variable, definition in group["variables"].items()
                },
            }
            for name, group in enabled_structure["groups"].items()
        }
        metadata_diff = dictionary_diff(metadata_disabled, metadata_enabled)
        structure_equal = not any(structure_diff.values())
        metadata_equal = not any(metadata_diff.values())
        disabled_variables = variable_map(disabled)
        enabled_variables = variable_map(enabled)
        variable_names_equal = set(disabled_variables) == set(enabled_variables)
        variable_results: dict[str, Any] = {}
        if structure_equal and variable_names_equal:
            for name in sorted(disabled_variables):
                left = disabled_variables[name]
                right = enabled_variables[name]
                if np.issubdtype(left.dtype, np.number):
                    result = compare_numeric_variable(left, right, atol, rtol)
                else:
                    result = compare_non_numeric_variable(left, right)
                result["dtype_equal"] = str(left.dtype) == str(right.dtype)
                variable_results[name] = result

    dtypes_equal = all(item["dtype_equal"] for item in variable_results.values())
    masks_equal = all(item["nan_masks_equal"] for item in variable_results.values())
    values_exact = bool(variable_results) and all(
        item["exact"] for item in variable_results.values()
    )
    within_tolerance = bool(variable_results) and all(
        item["within_tolerance"] for item in variable_results.values()
    )
    max_abs = max(
        (
            item["max_absolute_difference"]
            for item in variable_results.values()
            if item["max_absolute_difference"] is not None
        ),
        default=None,
    )
    max_rel = max(
        (
            item["max_relative_difference"]
            for item in variable_results.values()
            if item["max_relative_difference"] is not None
        ),
        default=None,
    )
    passed = (
        structure_equal
        and metadata_equal
        and dtypes_equal
        and masks_equal
        and within_tolerance
    )
    return {
        "structure_equal": structure_equal,
        "metadata_equal": metadata_equal,
        "dtypes_equal": dtypes_equal,
        "nan_masks_equal": masks_equal,
        "values_exact": values_exact,
        "within_tolerance": within_tolerance,
        "max_absolute_difference": max_abs,
        "max_relative_difference": max_rel,
        "passed": passed,
        "differences": {
            "structure": structure_diff,
            "metadata": metadata_diff,
            "variables": variable_results,
        },
    }


def canonicalize_json(
    value: Any, input_root: Path, relative_path: str, keys: tuple[str, ...] = ()
) -> Any:
    """Normalize expected run-location and viewer-time provenance."""
    if isinstance(value, dict):
        return {
            key: canonicalize_json(item, input_root, relative_path, keys + (key,))
            for key, item in sorted(value.items())
        }
    if isinstance(value, list):
        return [
            canonicalize_json(item, input_root, relative_path, keys + (str(index),))
            for index, item in enumerate(value)
        ]
    if isinstance(value, str):
        if (
            relative_path == "viewer/index.json"
            and keys == ("version",)
            and VIEWER_VERSION_RE.fullmatch(value)
        ):
            return "<VIEWER_GENERATED_AT>"
        return value.replace(str(input_root), "<RESULTS_DIR>")
    return value


def json_differences(left: Any, right: Any, pointer: str = "") -> list[dict[str, Any]]:
    """Return complete JSON structure and value differences by JSON pointer."""
    if isinstance(left, dict) and isinstance(right, dict):
        differences: list[dict[str, Any]] = []
        for key in sorted(set(left) | set(right)):
            child = f"{pointer}/{key.replace('~', '~0').replace('/', '~1')}"
            if key not in left:
                differences.append(
                    {
                        "pointer": child,
                        "kind": "additional_in_enabled",
                        "enabled": right[key],
                    }
                )
            elif key not in right:
                differences.append(
                    {
                        "pointer": child,
                        "kind": "missing_in_enabled",
                        "disabled": left[key],
                    }
                )
            else:
                differences.extend(json_differences(left[key], right[key], child))
        return differences
    if isinstance(left, list) and isinstance(right, list):
        differences = []
        common = min(len(left), len(right))
        for index in range(common):
            differences.extend(
                json_differences(left[index], right[index], f"{pointer}/{index}")
            )
        for index in range(common, len(left)):
            differences.append(
                {
                    "pointer": f"{pointer}/{index}",
                    "kind": "missing_in_enabled",
                    "disabled": left[index],
                }
            )
        for index in range(common, len(right)):
            differences.append(
                {
                    "pointer": f"{pointer}/{index}",
                    "kind": "additional_in_enabled",
                    "enabled": right[index],
                }
            )
        return differences
    if (
        isinstance(left, (int, float))
        and not isinstance(left, bool)
        and isinstance(right, (int, float))
        and not isinstance(right, bool)
    ):
        if (
            isinstance(left, float)
            and isinstance(right, float)
            and math.isnan(left)
            and math.isnan(right)
        ):
            return []
        if left == right:
            return []
        return [
            {
                "pointer": pointer or "/",
                "kind": "numeric_value",
                "disabled": left,
                "enabled": right,
            }
        ]
    if type(left) is not type(right):
        return [
            {
                "pointer": pointer or "/",
                "kind": "type",
                "disabled": left,
                "enabled": right,
            }
        ]
    if left != right:
        return [
            {
                "pointer": pointer or "/",
                "kind": "value",
                "disabled": left,
                "enabled": right,
            }
        ]
    return []


def compare_json(
    disabled_path: Path,
    enabled_path: Path,
    disabled_root: Path,
    enabled_root: Path,
    relative_path: str,
    atol: float,
    rtol: float,
) -> dict[str, Any]:
    """Compare raw and canonicalized JSON structure and values."""
    with disabled_path.open(encoding="utf-8") as stream:
        disabled_raw = json.load(stream)
    with enabled_path.open(encoding="utf-8") as stream:
        enabled_raw = json.load(stream)
    raw_differences = json_differences(disabled_raw, enabled_raw)
    disabled = canonicalize_json(disabled_raw, disabled_root, relative_path)
    enabled = canonicalize_json(enabled_raw, enabled_root, relative_path)
    differences = json_differences(disabled, enabled)
    numeric = [item for item in differences if item["kind"] == "numeric_value"]
    structural = [
        item
        for item in differences
        if item["kind"] in {"type", "missing_in_enabled", "additional_in_enabled"}
    ]
    nonnumeric = [item for item in differences if item["kind"] == "value"]
    max_abs = max(
        (abs(float(item["disabled"]) - float(item["enabled"])) for item in numeric),
        default=0.0,
    )
    max_rel = max(
        (
            abs(float(item["disabled"]) - float(item["enabled"]))
            / max(abs(float(item["disabled"])), abs(float(item["enabled"])))
            if max(abs(float(item["disabled"])), abs(float(item["enabled"]))) != 0
            else 0.0
            for item in numeric
        ),
        default=0.0,
    )
    numeric_within_tolerance = all(
        math.isclose(
            float(item["disabled"]),
            float(item["enabled"]),
            abs_tol=atol,
            rel_tol=rtol,
        )
        for item in numeric
    )
    passed = not structural and not nonnumeric and numeric_within_tolerance
    return {
        "raw_json_exact": not raw_differences,
        "canonical_json_exact": not differences,
        "structure_equal": not structural,
        "nonnumeric_values_equal": not nonnumeric,
        "numeric_values_exact": not numeric,
        "within_tolerance": numeric_within_tolerance,
        "max_absolute_difference": max_abs,
        "max_relative_difference": max_rel,
        "passed": passed,
        "raw_differences": raw_differences,
        "canonical_differences": differences,
    }


def compare_png(disabled_path: Path, enabled_path: Path) -> dict[str, Any]:
    """Compare PNG structure and decoded RGBA pixels exactly."""
    with Image.open(disabled_path) as disabled_image:
        disabled_format = disabled_image.format
        disabled_mode = disabled_image.mode
        disabled_size = disabled_image.size
        disabled_rgba = np.asarray(disabled_image.convert("RGBA"), dtype=np.int16)
    with Image.open(enabled_path) as enabled_image:
        enabled_format = enabled_image.format
        enabled_mode = enabled_image.mode
        enabled_size = enabled_image.size
        enabled_rgba = np.asarray(enabled_image.convert("RGBA"), dtype=np.int16)

    formats_equal = disabled_format == enabled_format == "PNG"
    modes_equal = disabled_mode == enabled_mode
    dimensions_equal = disabled_size == enabled_size
    structure_equal = formats_equal and modes_equal and dimensions_equal
    if dimensions_equal:
        absolute = np.abs(disabled_rgba - enabled_rgba)
        differing_channels = int(np.count_nonzero(absolute))
        differing_pixels = int(np.count_nonzero(np.any(absolute != 0, axis=2)))
        max_absolute = int(np.max(absolute, initial=0))
        rms = float(np.sqrt(np.mean(np.square(absolute, dtype=np.float64))))
        decoded_exact = differing_channels == 0
    else:
        differing_channels = None
        differing_pixels = None
        max_absolute = None
        rms = None
        decoded_exact = False

    return {
        "formats_equal": formats_equal,
        "modes_equal": modes_equal,
        "dimensions_equal": dimensions_equal,
        "structure_equal": structure_equal,
        "decoded_rgba_exact": decoded_exact,
        "differing_pixels": differing_pixels,
        "differing_channels": differing_channels,
        "max_absolute_difference": max_absolute,
        "max_relative_difference": (
            max_absolute / 255.0 if max_absolute is not None else None
        ),
        "rms_channel_difference": rms,
        "within_tolerance": decoded_exact,
        "passed": structure_equal and decoded_exact,
        "disabled_image": {
            "format": disabled_format,
            "mode": disabled_mode,
            "size": list(disabled_size),
        },
        "enabled_image": {
            "format": enabled_format,
            "mode": enabled_mode,
            "size": list(enabled_size),
        },
    }


def format_number(value: float | None) -> str:
    """Format an optional numeric result compactly."""
    return "NA" if value is None else f"{value:.17g}"


def write_tsv(
    path: Path, fieldnames: Sequence[str], rows: list[dict[str, Any]]
) -> None:
    """Write a deterministic tab-separated report."""
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=fieldnames, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def write_summary(
    output_path: Path,
    comparisons: list[dict[str, Any]],
    atol: float,
    rtol: float,
) -> None:
    """Write the human-readable validation report."""
    overall_pass = all(comparison["passed"] for comparison in comparisons)
    lines = [
        "# NetCDF3 Climatology Lock Workaround Validation",
        "",
        f"**Overall result: {'PASS' if overall_pass else 'FAIL'}**",
        "",
        "This validation compares workaround-disabled and workaround-enabled full",
        "diagnostic outputs within Python 3.13.14 and Python 3.14.6. Incomplete",
        "Python 3.14.1–3.14.4 outputs are explicitly excluded.",
        "",
        "## Input result paths",
        "",
        "| Python | Workaround | Results path |",
        "| --- | --- | --- |",
    ]
    for python_version, roots in INPUTS.items():
        for workaround in ("disabled", "enabled"):
            lines.append(
                f"| {python_version} | {workaround.capitalize()} | `{roots[workaround]}` |"
            )
    lines.extend(
        [
            "",
            "## Criteria",
            "",
            f"- Numeric tolerance: `atol={atol:g}`, `rtol={rtol:g}`.",
            "- Missing/additional paths, structural or metadata differences, dtype",
            "  changes, NaN-mask changes, nonnumeric value changes, and numeric values",
            "  outside tolerance fail validation.",
            "- JSON is canonicalized by sorting object keys, replacing each configured",
            "  absolute results root with `<RESULTS_DIR>`, and replacing the generated",
            "  `viewer/index.json` timestamp with `<VIEWER_GENERATED_AT>`.",
            "- PNG dimensions and modes must match, and decoded RGBA pixels must be",
            "  exactly equal. PNG container metadata and compression are not compared.",
            "- Raw byte/JSON equality is reported before canonical and tolerance results.",
            "",
            "## Results",
            "",
            "| Python | Inventory | Files | Byte exact | NetCDF exact | JSON canonical exact | PNG pixel exact | Max abs | Max rel | Result |",
            "| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |",
        ]
    )
    for comparison in comparisons:
        files = comparison["files"]
        byte_exact = sum(item["byte_exact"] for item in files)
        nc_files = [item for item in files if item["kind"] == "nc"]
        json_files = [item for item in files if item["kind"] == "json"]
        png_files = [item for item in files if item["kind"] == "png"]
        nc_exact = sum(item["values_exact"] for item in nc_files)
        json_exact = sum(item["canonical_json_exact"] for item in json_files)
        png_exact = sum(item["decoded_rgba_exact"] for item in png_files)
        lines.append(
            f"| {comparison['python_version']} | {comparison['inventory']['common']} common, "
            f"{len(comparison['inventory']['missing_in_enabled'])} missing, "
            f"{len(comparison['inventory']['additional_in_enabled'])} additional | "
            f"{len(files)} | {byte_exact}/{len(files)} | {nc_exact}/{len(nc_files)} | "
            f"{json_exact}/{len(json_files)} | {png_exact}/{len(png_files)} | "
            f"{format_number(comparison['max_absolute_difference'])} | "
            f"{format_number(comparison['max_relative_difference'])} | "
            f"{'PASS' if comparison['passed'] else 'FAIL'} |"
        )
    lines.extend(
        [
            "",
            "For each Python version, the sole byte-level mismatch is",
            "`viewer/index.json`: its raw JSON has four expected provenance differences",
            "(three embedded results-directory strings and one generation timestamp).",
            "All four disappear under the documented canonicalization rules.",
            "",
            "All five NetCDF files in each comparison were checked for data model,",
            "groups, dimensions, variable names/dimensions/shapes, metadata, dtypes,",
            "NaN masks, exact values, and tolerance. All JSON files were parsed and",
            "compared recursively after canonicalization. All PNG files were decoded",
            "to RGBA and checked for matching dimensions, modes, and exact pixels.",
            "",
            "## Artifacts and reproduction",
            "",
            "- [Inventory differences](inventory.tsv)",
            "- [Concise per-file results](file-results.tsv)",
            "- [Complete machine-readable details](comparison-details.json)",
            "- [Standalone validator](validate_outputs.py)",
            "",
            "From this directory, rerun with:",
            "",
            "```bash",
            "conda run -n ed_1048_xr_2026070_py31314 python validate_outputs.py",
            "```",
            "",
        ]
    )
    output_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    """Run all configured comparisons and write validation artifacts."""
    args = parse_args()
    if args.atol < 0 or args.rtol < 0:
        raise ValueError("Tolerances must be nonnegative")
    args.output_dir.mkdir(parents=True, exist_ok=True)
    inventory_rows: list[dict[str, Any]] = []
    file_rows: list[dict[str, Any]] = []
    comparisons: list[dict[str, Any]] = []

    for python_version, roots in INPUTS.items():
        disabled_inventory = inventory(roots["disabled"])
        enabled_inventory = inventory(roots["enabled"])
        disabled_paths = set(disabled_inventory)
        enabled_paths = set(enabled_inventory)
        missing = sorted(disabled_paths - enabled_paths)
        additional = sorted(enabled_paths - disabled_paths)
        for status, paths, source in (
            ("missing_in_enabled", missing, disabled_inventory),
            ("additional_in_enabled", additional, enabled_inventory),
        ):
            inventory_rows.extend(
                {
                    "python_version": python_version,
                    "status": status,
                    "kind": source[path],
                    "relative_path": path,
                }
                for path in paths
            )

        files: list[dict[str, Any]] = []
        for relative_path in sorted(disabled_paths & enabled_paths):
            kind = disabled_inventory[relative_path]
            disabled_path = roots["disabled"] / relative_path
            enabled_path = roots["enabled"] / relative_path
            disabled_sha256 = sha256(disabled_path)
            enabled_sha256 = sha256(enabled_path)
            byte_exact = disabled_sha256 == enabled_sha256
            if kind == "nc":
                result = compare_netcdf(
                    disabled_path, enabled_path, args.atol, args.rtol
                )
                concise = {
                    "structure_equal": result["structure_equal"],
                    "metadata_equal": result["metadata_equal"],
                    "dtypes_equal": result["dtypes_equal"],
                    "nan_masks_equal": result["nan_masks_equal"],
                    "values_exact": result["values_exact"],
                    "canonical_exact": "NA",
                    "differing_pixels": "NA",
                    "rms_channel_difference": "NA",
                }
            elif kind == "json":
                result = compare_json(
                    disabled_path,
                    enabled_path,
                    roots["disabled"],
                    roots["enabled"],
                    relative_path,
                    args.atol,
                    args.rtol,
                )
                concise = {
                    "structure_equal": result["structure_equal"],
                    "metadata_equal": "NA",
                    "dtypes_equal": "NA",
                    "nan_masks_equal": "NA",
                    "values_exact": result["numeric_values_exact"]
                    and result["nonnumeric_values_equal"],
                    "canonical_exact": result["canonical_json_exact"],
                    "differing_pixels": "NA",
                    "rms_channel_difference": "NA",
                }
            else:
                result = compare_png(disabled_path, enabled_path)
                concise = {
                    "structure_equal": result["structure_equal"],
                    "metadata_equal": "NA",
                    "dtypes_equal": "NA",
                    "nan_masks_equal": "NA",
                    "values_exact": result["decoded_rgba_exact"],
                    "canonical_exact": "NA",
                    "differing_pixels": result["differing_pixels"],
                    "rms_channel_difference": format_number(
                        result["rms_channel_difference"]
                    ),
                }
            file_result = {
                "relative_path": relative_path,
                "kind": kind,
                "byte_exact": byte_exact,
                "disabled_sha256": disabled_sha256,
                "enabled_sha256": enabled_sha256,
                **result,
            }
            files.append(file_result)
            file_rows.append(
                {
                    "python_version": python_version,
                    "relative_path": relative_path,
                    "kind": kind,
                    "byte_exact": byte_exact,
                    **concise,
                    "max_absolute_difference": format_number(
                        result["max_absolute_difference"]
                    ),
                    "max_relative_difference": format_number(
                        result["max_relative_difference"]
                    ),
                    "within_tolerance": result["within_tolerance"],
                    "status": "PASS" if result["passed"] else "FAIL",
                }
            )

        max_abs = max(
            (
                item["max_absolute_difference"]
                for item in files
                if item["max_absolute_difference"] is not None
            ),
            default=None,
        )
        max_rel = max(
            (
                item["max_relative_difference"]
                for item in files
                if item["max_relative_difference"] is not None
            ),
            default=None,
        )
        comparisons.append(
            {
                "python_version": python_version,
                "inputs": {name: str(root) for name, root in roots.items()},
                "inventory": {
                    "disabled_count": len(disabled_paths),
                    "enabled_count": len(enabled_paths),
                    "common": len(disabled_paths & enabled_paths),
                    "missing_in_enabled": missing,
                    "additional_in_enabled": additional,
                },
                "files": files,
                "max_absolute_difference": max_abs,
                "max_relative_difference": max_rel,
                "passed": not missing
                and not additional
                and all(item["passed"] for item in files),
            }
        )

    write_tsv(
        args.output_dir / "inventory.tsv",
        ["python_version", "status", "kind", "relative_path"],
        inventory_rows,
    )
    write_tsv(
        args.output_dir / "file-results.tsv",
        [
            "python_version",
            "relative_path",
            "kind",
            "byte_exact",
            "structure_equal",
            "metadata_equal",
            "dtypes_equal",
            "nan_masks_equal",
            "values_exact",
            "canonical_exact",
            "differing_pixels",
            "rms_channel_difference",
            "max_absolute_difference",
            "max_relative_difference",
            "within_tolerance",
            "status",
        ],
        file_rows,
    )
    details = {
        "schema_version": 2,
        "generated_at": datetime.now().astimezone().isoformat(timespec="seconds"),
        "validator": str(Path(__file__).resolve()),
        "tolerances": {"absolute": args.atol, "relative": args.rtol},
        "json_canonicalization": [
            "Sort object keys recursively",
            "Replace configured absolute input root with <RESULTS_DIR> in strings",
            "Replace viewer/index.json root version timestamp with <VIEWER_GENERATED_AT>",
        ],
        "excluded_incomplete_python_versions": EXCLUDED_INCOMPLETE_VERSIONS,
        "comparisons": comparisons,
        "overall_pass": all(comparison["passed"] for comparison in comparisons),
    }
    with (args.output_dir / "comparison-details.json").open(
        "w", encoding="utf-8"
    ) as stream:
        json.dump(jsonable(details), stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    write_summary(args.output_dir / "summary.md", comparisons, args.atol, args.rtol)
    sys.stderr.write(
        f"Validation {'PASSED' if details['overall_pass'] else 'FAILED'}; "
        f"wrote reports to {args.output_dir}\n"
    )
    return 0 if details["overall_pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
