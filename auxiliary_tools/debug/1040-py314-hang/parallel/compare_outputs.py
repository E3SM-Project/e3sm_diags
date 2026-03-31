#!/usr/bin/env python3
"""Compare PNG outputs across three public HTML result trees.

This script walks the configured output directories recursively, counts all PNG
files, compares matching images pairwise with Pillow, and writes a CSV with one
row per relative PNG path.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from PIL import Image, ImageChops, ImageStat
from tqdm import tqdm

from e3sm_diags.logger import _setup_child_logger


logger = _setup_child_logger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
OUTPUT_CSV = SCRIPT_DIR / "compare_outputs_metrics.csv"


@dataclass(frozen=True)
class OutputRoot:
    """Metadata for a PNG output tree."""

    label: str
    path: Path


ROOTS: tuple[OutputRoot, ...] = (
    OutputRoot(
        label="py313_main",
        path=Path("/lcrc/group/e3sm/public_html/ac.tvo/1040-py314-hang-tom-py313-main"),
    ),
    OutputRoot(
        label="py313",
        path=Path("/lcrc/group/e3sm/public_html/ac.tvo/1040-py314-hang-tom-py313"),
    ),
    OutputRoot(
        label="py314",
        path=Path("/lcrc/group/e3sm/public_html/ac.tvo/1040-py314-hang-tom-py314"),
    ),
)


def list_pngs(root: Path) -> dict[str, Path]:
    """Return a mapping of relative PNG path to absolute file path."""
    pngs: dict[str, Path] = {}

    for path in root.rglob("*"):
        if path.is_file() and path.suffix.lower() == ".png":
            pngs[path.relative_to(root).as_posix()] = path

    return pngs


def _blank_metrics(status: str) -> dict[str, str]:
    return {
        "status": status,
        "compared": "False",
        "identical": "",
        "left_mode": "",
        "right_mode": "",
        "left_size": "",
        "right_size": "",
        "pixel_count": "",
        "pixel_diff_count": "",
        "pixel_diff_fraction": "",
        "mean_abs_diff": "",
        "rms_diff": "",
    }


def _count_diff_pixels(diff_image: Image.Image) -> int:
    bbox = diff_image.getbbox()
    if bbox is None:
        return 0

    return sum(
        diff_image.crop(bbox)
        .point(lambda value: 255 if value else 0)
        .convert("L")
        .point(bool)
        .getdata()
    )


def compare_pngs(left_path: Path | None, right_path: Path | None) -> dict[str, str]:
    """Compare two PNGs and return CSV-ready metrics."""
    if left_path is None and right_path is None:
        return _blank_metrics("missing_both")
    if left_path is None:
        return _blank_metrics("missing_left")
    if right_path is None:
        return _blank_metrics("missing_right")

    with Image.open(left_path) as left_raw, Image.open(right_path) as right_raw:
        left_mode = left_raw.mode
        right_mode = right_raw.mode
        left_size = left_raw.size
        right_size = right_raw.size

        if left_size != right_size:
            metrics = _blank_metrics("size_mismatch")
            metrics["left_mode"] = left_mode
            metrics["right_mode"] = right_mode
            metrics["left_size"] = f"{left_size[0]}x{left_size[1]}"
            metrics["right_size"] = f"{right_size[0]}x{right_size[1]}"
            return metrics

        left_image = left_raw.convert("RGB")
        right_image = right_raw.convert("RGB")

    diff = ImageChops.difference(left_image, right_image)
    diff_pixels = _count_diff_pixels(diff)
    pixel_count = left_image.size[0] * left_image.size[1]
    stat = ImageStat.Stat(diff)
    mean_abs_diff = sum(stat.mean) / len(stat.mean)
    rms_diff = sum(stat.rms) / len(stat.rms)

    return {
        "status": "compared",
        "compared": "True",
        "identical": str(diff.getbbox() is None),
        "left_mode": left_mode,
        "right_mode": right_mode,
        "left_size": f"{left_size[0]}x{left_size[1]}",
        "right_size": f"{right_size[0]}x{right_size[1]}",
        "pixel_count": str(pixel_count),
        "pixel_diff_count": str(diff_pixels),
        "pixel_diff_fraction": f"{diff_pixels / pixel_count:.12f}",
        "mean_abs_diff": f"{mean_abs_diff:.12f}",
        "rms_diff": f"{rms_diff:.12f}",
    }


def build_fieldnames() -> list[str]:
    """Build the CSV header."""
    fieldnames = ["relative_path"]

    for root in ROOTS:
        fieldnames.extend(
            [
                f"exists_{root.label}",
                f"png_count_{root.label}",
                f"file_bytes_{root.label}",
                f"path_{root.label}",
            ]
        )

    for left_root, right_root in combinations(ROOTS, 2):
        prefix = f"{left_root.label}_vs_{right_root.label}"
        fieldnames.extend(
            [
                f"{prefix}_status",
                f"{prefix}_compared",
                f"{prefix}_identical",
                f"{prefix}_left_mode",
                f"{prefix}_right_mode",
                f"{prefix}_left_size",
                f"{prefix}_right_size",
                f"{prefix}_pixel_count",
                f"{prefix}_pixel_diff_count",
                f"{prefix}_pixel_diff_fraction",
                f"{prefix}_mean_abs_diff",
                f"{prefix}_rms_diff",
            ]
        )

    return fieldnames


def main() -> int:
    """Run the PNG inventory and comparison workflow."""
    logger.info("Starting PNG comparison workflow")
    logger.info("Metrics CSV will be written to %s", OUTPUT_CSV)

    for root in ROOTS:
        logger.info("Scanning PNGs under %s: %s", root.label, root.path)

    inventory = {root.label: list_pngs(root.path) for root in ROOTS}
    counts = {root.label: len(inventory[root.label]) for root in ROOTS}
    relative_paths = sorted(
        {
            rel_path
            for pngs_by_path in inventory.values()
            for rel_path in pngs_by_path
        }
    )

    logger.info("PNG counts by root:")
    for root in ROOTS:
        logger.info("  %s: %s", root.label, counts[root.label])

    logger.info(
        "Found %s unique relative PNG paths across %s roots",
        len(relative_paths),
        len(ROOTS),
    )
    logger.info(
        "Running %s pairwise comparisons per PNG: %s",
        3,
        ", ".join(
            f"{left_root.label}_vs_{right_root.label}"
            for left_root, right_root in combinations(ROOTS, 2)
        ),
    )
    logger.info("Writing CSV rows")

    with OUTPUT_CSV.open("w", newline="", encoding="utf-8") as csv_file:
        writer = csv.DictWriter(csv_file, fieldnames=build_fieldnames())
        writer.writeheader()

        for relative_path in tqdm(
            relative_paths,
            desc="Comparing PNGs",
            unit="image",
        ):
            row: dict[str, str] = {"relative_path": relative_path}

            for root in ROOTS:
                path = inventory[root.label].get(relative_path)
                row[f"exists_{root.label}"] = str(path is not None)
                row[f"png_count_{root.label}"] = str(counts[root.label])
                row[f"file_bytes_{root.label}"] = str(path.stat().st_size) if path else ""
                row[f"path_{root.label}"] = str(path) if path else ""

            for left_root, right_root in combinations(ROOTS, 2):
                prefix = f"{left_root.label}_vs_{right_root.label}"
                metrics = compare_pngs(
                    inventory[left_root.label].get(relative_path),
                    inventory[right_root.label].get(relative_path),
                )
                for key, value in metrics.items():
                    row[f"{prefix}_{key}"] = value

            writer.writerow(row)

    logger.info("Wrote comparison metrics CSV to %s", OUTPUT_CSV)
    logger.info("Compared %s unique relative PNG paths", len(relative_paths))
    logger.info("PNG comparison workflow complete")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
