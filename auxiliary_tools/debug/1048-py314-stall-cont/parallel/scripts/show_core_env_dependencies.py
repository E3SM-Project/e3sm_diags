#!/usr/bin/env python3
"""Show core package versions for the Python 3.14 debug environments."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

from e3sm_diags.logger import _setup_child_logger


logger = _setup_child_logger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_CSV = SCRIPT_DIR / "core_environment_dependencies.csv"
ENV_NAMES = (
    "ed_main_py313",
    "ed_main_py314",
    "ed_main_py314_xr2025101",
    "ed_main_py314_xr2025110",
    "ed_main_py314_xr2025120",
    "ed_main_py314_xr2026010",
    "ed_1048_py314",
)
CONDA_EXE = Path("/home/ac.tvo/miniforge3/bin/conda")
CORE_PACKAGES = (
    ("python", ("python",)),
    ("xarray", ("xarray",)),
    ("uxarray", ("uxarray",)),
    ("numpy", ("numpy",)),
    ("pandas", ("pandas",)),
    ("dask", ("dask",)),
    ("matplotlib", ("matplotlib", "matplotlib-base")),
    ("cartopy", ("cartopy",)),
)


def load_packages(env_name: str) -> dict[str, dict[str, str]]:
    """Load package metadata from ``conda list --json``."""
    command = [str(CONDA_EXE), "list", "-n", env_name, "--json"]
    proc = subprocess.run(command, check=True, capture_output=True, text=True)
    raw_packages = json.loads(proc.stdout)

    return {
        pkg["name"]: {
            "version": pkg.get("version", ""),
            "channel": pkg.get("channel", ""),
        }
        for pkg in raw_packages
    }


def resolve_package(
    env_packages: dict[str, dict[str, str]],
    candidates: tuple[str, ...],
) -> dict[str, str]:
    """Resolve the first matching package from the candidate names."""
    for package_name in candidates:
        package = env_packages.get(package_name)
        if package is not None:
            return {
                "resolved_package": package_name,
                "version": package.get("version", ""),
                "channel": package.get("channel", ""),
            }

    return {
        "resolved_package": "",
        "version": "",
        "channel": "",
    }


def build_rows(
    packages_by_env: dict[str, dict[str, dict[str, str]]]
) -> list[dict[str, str]]:
    """Build CSV rows for the configured core packages."""
    rows: list[dict[str, str]] = []

    for display_name, candidates in CORE_PACKAGES:
        row = {"package": display_name}

        for env_name in ENV_NAMES:
            resolved = resolve_package(packages_by_env[env_name], candidates)
            row[f"{env_name}_resolved_package"] = resolved["resolved_package"]
            row[f"{env_name}_version"] = resolved["version"]
            row[f"{env_name}_channel"] = resolved["channel"]

        rows.append(row)

    return rows


def write_rows(output_csv: Path, rows: list[dict[str, str]]) -> None:
    """Write the core package summary to CSV."""
    fieldnames = ["package"]
    for env_name in ENV_NAMES:
        fieldnames.extend(
            [
                f"{env_name}_resolved_package",
                f"{env_name}_version",
                f"{env_name}_channel",
            ]
        )

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def format_summary_value(row: dict[str, str], env_name: str) -> str:
    """Format a human-readable summary value for one environment."""
    resolved_package = row[f"{env_name}_resolved_package"]
    version = row[f"{env_name}_version"]
    channel = row[f"{env_name}_channel"]

    if not resolved_package:
        return "MISSING"

    summary = f"{version} [{channel or 'unknown'}]"
    if resolved_package != row["package"]:
        summary = f"{summary} via {resolved_package}"

    return summary


def print_summary(rows: list[dict[str, str]]) -> None:
    """Print a concise summary to stdout."""
    print("Core dependency summary")
    for row in rows:
        print(row["package"])
        for env_name in ENV_NAMES:
            print(f"  {env_name}: {format_summary_value(row, env_name)}")


def main() -> int:
    """Load, write, and print the core dependency summary."""
    logger.info("Loading conda package lists for %s", ", ".join(ENV_NAMES))
    packages_by_env = {
        env_name: load_packages(env_name)
        for env_name in ENV_NAMES
    }
    rows = build_rows(packages_by_env)
    logger.info("Writing core environment dependencies to %s", OUTPUT_CSV)
    write_rows(OUTPUT_CSV, rows)
    print_summary(rows)
    logger.info("Core environment dependency summary complete")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
