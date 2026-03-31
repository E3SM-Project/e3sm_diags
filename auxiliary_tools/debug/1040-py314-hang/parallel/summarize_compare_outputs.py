#!/usr/bin/env python3
"""Summarize aggregate PNG comparison metrics from compare_outputs_metrics.csv."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

from e3sm_diags.logger import _setup_child_logger


logger = _setup_child_logger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_CSV = SCRIPT_DIR / "compare_outputs_metrics.csv"
DEFAULT_OUTPUT_CSV = SCRIPT_DIR / "compare_outputs_summary.csv"
DEFAULT_DIFFS_CSV = SCRIPT_DIR / "compare_outputs_differences.csv"
ROOT_LABELS = ("py313_main", "py313", "py314")
PAIR_LABELS = (
    "py313_main_vs_py313",
    "py313_main_vs_py314",
    "py313_vs_py314",
)


def _to_bool(value: str) -> bool:
    return value == "True"


def _to_int(value: str) -> int:
    if value == "":
        return 0
    return int(value)


def _to_float(value: str) -> float:
    if value == "":
        return 0.0
    return float(value)


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Read compare_outputs_metrics.csv and compute aggregate PNG "
            "comparison metrics."
        )
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DEFAULT_INPUT_CSV,
        help="Path to the detailed comparison metrics CSV",
    )
    parser.add_argument(
        "--output-csv",
        type=Path,
        default=DEFAULT_OUTPUT_CSV,
        help="Path to write the aggregate summary CSV",
    )
    parser.add_argument(
        "--diffs-csv",
        type=Path,
        default=DEFAULT_DIFFS_CSV,
        help="Path to write the per-file difference CSV",
    )

    return parser.parse_args()


def init_summary_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Build aggregate summary rows from detailed comparison rows."""
    summary_rows: list[dict[str, str]] = []

    summary_rows.append(
        {
            "category": "overall",
            "name": "unique_relative_png_paths",
            "value": str(len(rows)),
        }
    )

    for root_label in ROOT_LABELS:
        exists_key = f"exists_{root_label}"
        total_present = sum(_to_bool(row[exists_key]) for row in rows)
        summary_rows.extend(
            [
                {
                    "category": "branch",
                    "name": f"{root_label}_present_files",
                    "value": str(total_present),
                },
                {
                    "category": "branch",
                    "name": f"{root_label}_missing_files",
                    "value": str(len(rows) - total_present),
                },
            ]
        )

    for root_label in ROOT_LABELS[1:]:
        summary_rows.extend(
            [
                {
                    "category": "vs_main",
                    "name": f"{root_label}_missing_compared_to_main",
                    "value": str(
                        sum(
                            _to_bool(row["exists_py313_main"])
                            and not _to_bool(row[f"exists_{root_label}"])
                            for row in rows
                        )
                    ),
                },
                {
                    "category": "vs_main",
                    "name": f"{root_label}_extra_compared_to_main",
                    "value": str(
                        sum(
                            not _to_bool(row["exists_py313_main"])
                            and _to_bool(row[f"exists_{root_label}"])
                            for row in rows
                        )
                    ),
                },
            ]
        )

    for pair_label in PAIR_LABELS:
        status_key = f"{pair_label}_status"
        compared_key = f"{pair_label}_compared"
        identical_key = f"{pair_label}_identical"
        pixel_diff_key = f"{pair_label}_pixel_diff_count"
        pixel_fraction_key = f"{pair_label}_pixel_diff_fraction"

        compared_rows = [row for row in rows if _to_bool(row[compared_key])]
        different_rows = [
            row for row in compared_rows if not _to_bool(row[identical_key])
        ]

        summary_rows.extend(
            [
                {
                    "category": "pair",
                    "name": f"{pair_label}_compared_files",
                    "value": str(len(compared_rows)),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_identical_files",
                    "value": str(
                        sum(_to_bool(row[identical_key]) for row in compared_rows)
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_different_files",
                    "value": str(len(different_rows)),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_missing_left",
                    "value": str(
                        sum(row[status_key] == "missing_left" for row in rows)
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_missing_right",
                    "value": str(
                        sum(row[status_key] == "missing_right" for row in rows)
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_size_mismatch",
                    "value": str(
                        sum(row[status_key] == "size_mismatch" for row in rows)
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_total_pixel_diff_count",
                    "value": str(
                        sum(_to_int(row[pixel_diff_key]) for row in compared_rows)
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_max_pixel_diff_fraction",
                    "value": "{:.12f}".format(
                        max(
                            (_to_float(row[pixel_fraction_key]) for row in compared_rows),
                            default=0.0,
                        )
                    ),
                },
                {
                    "category": "pair",
                    "name": f"{pair_label}_mean_pixel_diff_fraction",
                    "value": "{:.12f}".format(
                        (
                            sum(
                                _to_float(row[pixel_fraction_key])
                                for row in compared_rows
                            )
                            / len(compared_rows)
                        )
                        if compared_rows
                        else 0.0
                    ),
                },
            ]
        )

    return summary_rows


def read_rows(input_csv: Path) -> list[dict[str, str]]:
    """Read the detailed comparison CSV."""
    with input_csv.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def build_diff_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    """Build a focused per-file CSV for mismatches and missing files."""
    diff_rows: list[dict[str, str]] = []

    for row in rows:
        for pair_label in PAIR_LABELS:
            status = row[f"{pair_label}_status"]
            compared = _to_bool(row[f"{pair_label}_compared"])
            identical = row[f"{pair_label}_identical"]

            include_row = status != "compared" or (
                compared and identical == "False"
            )
            if not include_row:
                continue

            left_label, right_label = pair_label.split("_vs_")
            diff_rows.append(
                {
                    "pair": pair_label,
                    "relative_path": row["relative_path"],
                    "status": status,
                    "compared": row[f"{pair_label}_compared"],
                    "identical": identical,
                    "pixel_diff_count": row[f"{pair_label}_pixel_diff_count"],
                    "pixel_diff_fraction": row[f"{pair_label}_pixel_diff_fraction"],
                    "mean_abs_diff": row[f"{pair_label}_mean_abs_diff"],
                    "rms_diff": row[f"{pair_label}_rms_diff"],
                    "left_size": row[f"{pair_label}_left_size"],
                    "right_size": row[f"{pair_label}_right_size"],
                    "left_path": row[f"path_{left_label}"],
                    "right_path": row[f"path_{right_label}"],
                }
            )

    diff_rows.sort(
        key=lambda item: (
            item["pair"],
            -_to_float(item["pixel_diff_fraction"]),
            item["relative_path"],
        )
    )

    return diff_rows


def write_summary(output_csv: Path, summary_rows: list[dict[str, str]]) -> None:
    """Write the aggregate summary CSV."""
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["category", "name", "value"])
        writer.writeheader()
        writer.writerows(summary_rows)


def write_diffs(output_csv: Path, diff_rows: list[dict[str, str]]) -> None:
    """Write the per-file difference CSV."""
    with output_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pair",
                "relative_path",
                "status",
                "compared",
                "identical",
                "pixel_diff_count",
                "pixel_diff_fraction",
                "mean_abs_diff",
                "rms_diff",
                "left_size",
                "right_size",
                "left_path",
                "right_path",
            ],
        )
        writer.writeheader()
        writer.writerows(diff_rows)


def print_summary(summary_rows: list[dict[str, str]]) -> None:
    """Print a concise human-readable summary."""
    summary_map = {row["name"]: row["value"] for row in summary_rows}

    print("Branch totals:")
    for root_label in ROOT_LABELS:
        print(
            "  {}: present={}, missing={}".format(
                root_label,
                summary_map[f"{root_label}_present_files"],
                summary_map[f"{root_label}_missing_files"],
            )
        )

    print("\nCompared to main:")
    for root_label in ROOT_LABELS[1:]:
        print(
            "  {}: missing_vs_main={}, extra_vs_main={}".format(
                root_label,
                summary_map[f"{root_label}_missing_compared_to_main"],
                summary_map[f"{root_label}_extra_compared_to_main"],
            )
        )

    print("\nPairwise image diffs:")
    for pair_label in PAIR_LABELS:
        print(
            "  {}: compared={}, identical={}, different={}, total_pixel_diff={}, "
            "max_fraction={}".format(
                pair_label,
                summary_map[f"{pair_label}_compared_files"],
                summary_map[f"{pair_label}_identical_files"],
                summary_map[f"{pair_label}_different_files"],
                summary_map[f"{pair_label}_total_pixel_diff_count"],
                summary_map[f"{pair_label}_max_pixel_diff_fraction"],
            )
        )


def main() -> int:
    """Load the detailed metrics CSV and produce aggregate summaries."""
    args = parse_args()

    logger.info("Reading detailed metrics from %s", args.input_csv)
    rows = read_rows(args.input_csv)
    logger.info("Loaded %s detailed rows", len(rows))

    summary_rows = init_summary_rows(rows)
    diff_rows = build_diff_rows(rows)
    logger.info("Writing aggregate summary to %s", args.output_csv)
    write_summary(args.output_csv, summary_rows)
    logger.info("Writing per-file differences to %s", args.diffs_csv)
    write_diffs(args.diffs_csv, diff_rows)
    print_summary(summary_rows)
    print("\nDetailed diff rows written: {}".format(len(diff_rows)))
    logger.info("Aggregate summary workflow complete")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
