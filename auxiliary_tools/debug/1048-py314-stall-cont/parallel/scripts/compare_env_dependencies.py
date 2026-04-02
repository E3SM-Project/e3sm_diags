#!/usr/bin/env python3
"""Compare conda package dependencies across the three debug environments."""

from __future__ import annotations

import csv
import json
import subprocess
from pathlib import Path

from e3sm_diags.logger import _setup_child_logger


logger = _setup_child_logger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_CSV = SCRIPT_DIR / "environment_dependency_differences.csv"
ENV_NAMES = (
    "ed_main_py314",
    "ed_1048_py314",
)
CONDA_EXE = Path("/home/ac.tvo/miniforge3/bin/conda")


def load_packages(env_name: str) -> dict[str, dict[str, str]]:
    """Load package metadata from `conda list --json`."""
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


def build_rows(
    packages_by_env: dict[str, dict[str, dict[str, str]]]
) -> list[dict[str, str]]:
    """Build CSV rows for packages that differ across environments."""
    all_packages = sorted(
        {
            package_name
            for env_packages in packages_by_env.values()
            for package_name in env_packages
        }
    )
    rows: list[dict[str, str]] = []

    for package_name in all_packages:
        versions = []
        row = {"package": package_name}

        for env_name in ENV_NAMES:
            package = packages_by_env[env_name].get(package_name)
            version = package["version"] if package else ""
            channel = package["channel"] if package else ""
            row[f"{env_name}_version"] = version
            row[f"{env_name}_channel"] = channel
            versions.append((version, channel))

        if len(set(versions)) > 1:
            rows.append(row)

    return rows


def write_rows(output_csv: Path, rows: list[dict[str, str]]) -> None:
    """Write package differences to CSV."""
    fieldnames = ["package"]
    for env_name in ENV_NAMES:
        fieldnames.extend(
            [
                f"{env_name}_version",
                f"{env_name}_channel",
            ]
        )

    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def print_summary(rows: list[dict[str, str]]) -> None:
    """Print a concise summary to stdout."""
    print("Dependency differences: {}".format(len(rows)))
    for row in rows:
        print(row["package"])
        for env_name in ENV_NAMES:
            print(
                "  {}: {} [{}]".format(
                    env_name,
                    row[f"{env_name}_version"] or "MISSING",
                    row[f"{env_name}_channel"] or "MISSING",
                )
            )


def main() -> int:
    """Compare package sets across the configured environments."""
    logger.info("Loading conda package lists for %s", ", ".join(ENV_NAMES))
    packages_by_env = {
        env_name: load_packages(env_name)
        for env_name in ENV_NAMES
    }
    rows = build_rows(packages_by_env)
    logger.info("Writing environment dependency differences to %s", OUTPUT_CSV)
    write_rows(OUTPUT_CSV, rows)
    print_summary(rows)
    logger.info("Environment dependency comparison complete")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
