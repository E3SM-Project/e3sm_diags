#!/usr/bin/env python3
"""Render visual Pillow diff artifacts for rows in compare_outputs_differences.csv."""

from __future__ import annotations

import argparse
import csv
import shutil
from pathlib import Path

from PIL import Image, ImageChops, ImageDraw
from tqdm import tqdm

from e3sm_diags.logger import _setup_child_logger


logger = _setup_child_logger(__name__)

SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_INPUT_CSV = SCRIPT_DIR / "compare_outputs_differences.csv"
DEFAULT_OUTPUT_DIR = SCRIPT_DIR / "compare_output_diff_artifacts"
DEFAULT_MANIFEST_CSV = SCRIPT_DIR / "compare_output_diff_artifacts.csv"


def parse_args() -> argparse.Namespace:
    """Parse CLI arguments."""
    parser = argparse.ArgumentParser(
        description=(
            "Read compare_outputs_differences.csv and render Pillow diff "
            "artifacts for changed PNG pairs."
        )
    )
    parser.add_argument(
        "--input-csv",
        type=Path,
        default=DEFAULT_INPUT_CSV,
        help="Path to compare_outputs_differences.csv",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory where visual diff artifacts will be written",
    )
    parser.add_argument(
        "--manifest-csv",
        type=Path,
        default=DEFAULT_MANIFEST_CSV,
        help="CSV manifest describing the generated diff artifacts",
    )
    parser.add_argument(
        "--pair",
        default="",
        help="Optional pair filter such as py313_main_vs_py313",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=0,
        help="Optional maximum number of diff rows to render",
    )

    return parser.parse_args()


def read_rows(input_csv: Path) -> list[dict[str, str]]:
    """Read the differences CSV."""
    with input_csv.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def select_rows(rows: list[dict[str, str]], pair: str) -> list[dict[str, str]]:
    """Select only compared rows with actual image differences."""
    selected = [
        row
        for row in rows
        if row["status"] == "compared"
        and row["compared"] == "True"
        and row["identical"] == "False"
    ]
    if pair:
        selected = [row for row in selected if row["pair"] == pair]

    return selected


def build_overlay(base_image: Image.Image, diff_mask: Image.Image) -> Image.Image:
    """Create a red-highlight overlay that marks differing pixels."""
    overlay = base_image.convert("RGBA")
    mask = (
        diff_mask.convert("RGB")
        .point(lambda value: 255 if value else 0)
        .convert("L")
    )
    alpha_mask = mask.point(lambda value: 120 if value else 0)
    red_layer = Image.new("RGBA", overlay.size, (255, 0, 0, 0))
    red_layer.putalpha(alpha_mask)
    overlay.alpha_composite(red_layer)

    bbox = diff_mask.getbbox()
    if bbox is not None:
        draw = ImageDraw.Draw(overlay)
        draw.rectangle(bbox, outline=(255, 0, 0, 255), width=3)

    return overlay


def amplify_diff(diff_image: Image.Image, scale_factor: int = 16) -> Image.Image:
    """Amplify a diff image so tiny changes are easier to inspect."""
    return diff_image.point(
        lambda value: min(255, value * scale_factor)
    ).convert("RGB")


def render_row(output_dir: Path, row: dict[str, str]) -> dict[str, str]:
    """Render diff artifacts for a single CSV row."""
    relative_path = Path(row["relative_path"])
    pair_dir = output_dir / row["pair"] / relative_path.parent / relative_path.stem
    pair_dir.mkdir(parents=True, exist_ok=True)

    left_path = Path(row["left_path"])
    right_path = Path(row["right_path"])

    with Image.open(left_path) as left_raw, Image.open(right_path) as right_raw:
        left_image = left_raw.convert("RGB")
        right_image = right_raw.convert("RGB")
        diff_image = ImageChops.difference(left_image, right_image)

    overlay_image = build_overlay(right_image, diff_image)
    amplified_diff = amplify_diff(diff_image)
    bbox = diff_image.getbbox()
    bbox_text = ""
    if bbox is not None:
        bbox_text = "{},{},{},{}".format(*bbox)

    left_copy = pair_dir / "left.png"
    right_copy = pair_dir / "right.png"
    diff_raw_path = pair_dir / "diff_raw.png"
    diff_amplified_path = pair_dir / "diff_amplified.png"
    overlay_path = pair_dir / "overlay.png"

    shutil.copy2(left_path, left_copy)
    shutil.copy2(right_path, right_copy)
    diff_image.save(diff_raw_path, "PNG")
    amplified_diff.save(diff_amplified_path, "PNG")
    overlay_image.save(overlay_path, "PNG")

    return {
        "pair": row["pair"],
        "relative_path": row["relative_path"],
        "pixel_diff_count": row["pixel_diff_count"],
        "pixel_diff_fraction": row["pixel_diff_fraction"],
        "bbox": bbox_text,
        "left_copy": str(left_copy),
        "right_copy": str(right_copy),
        "diff_raw": str(diff_raw_path),
        "diff_amplified": str(diff_amplified_path),
        "overlay": str(overlay_path),
    }


def write_manifest(manifest_csv: Path, rows: list[dict[str, str]]) -> None:
    """Write a manifest of generated diff artifacts."""
    with manifest_csv.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "pair",
                "relative_path",
                "pixel_diff_count",
                "pixel_diff_fraction",
                "bbox",
                "left_copy",
                "right_copy",
                "diff_raw",
                "diff_amplified",
                "overlay",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    """Render visual diff artifacts for changed PNG pairs."""
    args = parse_args()

    logger.info("Reading diff rows from %s", args.input_csv)
    rows = read_rows(args.input_csv)
    selected_rows = select_rows(rows, args.pair)
    if args.limit > 0:
        selected_rows = selected_rows[: args.limit]

    logger.info("Rendering diff artifacts for %s rows", len(selected_rows))
    logger.info("Artifacts will be written under %s", args.output_dir)

    manifest_rows = []
    for row in tqdm(selected_rows, desc="Rendering diff artifacts", unit="diff"):
        manifest_rows.append(render_row(args.output_dir, row))

    logger.info("Writing artifact manifest to %s", args.manifest_csv)
    write_manifest(args.manifest_csv, manifest_rows)
    logger.info("Rendered %s diff artifact sets", len(manifest_rows))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
