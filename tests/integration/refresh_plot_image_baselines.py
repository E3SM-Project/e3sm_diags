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
    BASELINE_DIR,
    BASELINE_IMAGE_FILENAMES,
    render_lat_lon_plot_regression,
)


def refresh_lat_lon_plot_baselines(baseline_dir: str | Path) -> Path:
    baseline_path = Path(baseline_dir)
    baseline_path.mkdir(parents=True, exist_ok=True)

    with TemporaryDirectory() as temp_dir:
        _, generated_images = render_lat_lon_plot_regression(temp_dir)

        for source_path, filename in zip(generated_images, BASELINE_IMAGE_FILENAMES):
            shutil.copy2(source_path, baseline_path / filename)

    write_runtime_metadata(baseline_path / BASELINE_METADATA_FILENAME)

    return baseline_path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--baseline-dir",
        type=Path,
        default=BASELINE_DIR,
        help="Directory containing the committed targeted image-regression baselines.",
    )
    args = parser.parse_args()

    refresh_lat_lon_plot_baselines(args.baseline_dir)


if __name__ == "__main__":
    main()
