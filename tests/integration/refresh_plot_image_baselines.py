from __future__ import annotations

import argparse
import shutil
from pathlib import Path
from tempfile import TemporaryDirectory

from tests.integration.image_regression import (
    BASELINE_METADATA_FILENAME,
    write_runtime_metadata,
)
from tests.integration.plot_image_regression_case import (
    IMAGE_REGRESSION_CASES,
    IMAGE_REGRESSION_CASES_BY_ID,
    ImageRegressionCase,
)


def refresh_case_baselines(
    case: ImageRegressionCase, baseline_dir: str | Path | None = None
) -> Path:
    baseline_path = case.baseline_dir if baseline_dir is None else Path(baseline_dir)
    baseline_path.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory() as temp_dir:
        generated_images = case.render(temp_dir)

        for source_path, filename in zip(
            generated_images, case.expected_image_filenames, strict=True
        ):
            shutil.copy2(source_path, baseline_path / filename)

    metadata_path = (
        case.baseline_metadata_path
        if baseline_dir is None
        else baseline_path / BASELINE_METADATA_FILENAME
    )
    write_runtime_metadata(metadata_path)

    return baseline_path


def refresh_all_baselines() -> tuple[Path, ...]:
    return tuple(refresh_case_baselines(case) for case in IMAGE_REGRESSION_CASES)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--case",
        choices=tuple(IMAGE_REGRESSION_CASES_BY_ID),
        help="Refresh baselines for one targeted image-regression case.",
    )
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        help="Override baseline directory for single-case refreshes.",
    )
    args = parser.parse_args()

    if args.baseline_dir is not None and args.case is None:
        parser.error("--baseline-dir requires --case")

    if args.case is None:
        refresh_all_baselines()
        return

    refresh_case_baselines(
        IMAGE_REGRESSION_CASES_BY_ID[args.case],
        baseline_dir=args.baseline_dir,
    )


if __name__ == "__main__":
    main()
