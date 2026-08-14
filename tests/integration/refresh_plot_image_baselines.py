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
from tests.integration.utils import _compare_images


def refresh_case_baselines(
    case: ImageRegressionCase, baseline_dir: str | Path | None = None
) -> Path:
    baseline_path = case.baseline_dir if baseline_dir is None else Path(baseline_dir)
    baseline_path.mkdir(parents=True, exist_ok=True)
    has_changes = False

    with TemporaryDirectory() as temp_dir:
        generated_images = case.render(temp_dir)

        for source_path, filename in zip(
            generated_images, case.expected_image_filenames, strict=True
        ):
            expected_path = baseline_path / filename

            if not expected_path.exists():
                shutil.copy2(source_path, expected_path)
                has_changes = True
                continue

            mismatched_images: list[str] = []
            _compare_images(
                mismatched_images,
                filename,
                str(source_path),
                str(expected_path),
                diff_dir=str(Path(temp_dir) / "diff" / Path(filename).stem),
                mismatch_threshold=case.get_mismatch_threshold(filename),
            )

            if mismatched_images:
                shutil.copy2(source_path, expected_path)
                has_changes = True

    metadata_path = (
        case.baseline_metadata_path
        if baseline_dir is None
        else baseline_path / BASELINE_METADATA_FILENAME
    )
    if has_changes:
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
