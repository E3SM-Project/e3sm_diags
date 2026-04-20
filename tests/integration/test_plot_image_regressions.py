from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.integration.image_regression import (
    BASELINE_METADATA_FILENAME,
    assert_image_matches_baseline,
)
from tests.integration.plot_image_regression_case import (
    BASELINE_DIR,
    get_image_regression_artifact_dir,
    render_lat_lon_plot_regression,
)


@pytest.mark.image_regression
class TestLatLonPlotImageRegressions:
    def test_lat_lon_plot_matches_expected_baselines(self, tmp_path: Path):
        _, generated_images = render_lat_lon_plot_regression(tmp_path / "results")
        diff_artifact_dir = get_image_regression_artifact_dir(tmp_path)

        runtime_metadata_path = assert_image_matches_baseline(
            image_name="lat_lon_plot_regression.png",
            path_to_actual_png=generated_images[0],
            path_to_expected_png=BASELINE_DIR / "lat_lon_plot_regression.png",
            artifact_dir=diff_artifact_dir / "main_plot",
        )
        assert runtime_metadata_path.exists()

        assert_image_matches_baseline(
            image_name="lat_lon_plot_regression.2.png",
            path_to_actual_png=generated_images[1],
            path_to_expected_png=BASELINE_DIR / "lat_lon_plot_regression.2.png",
            artifact_dir=diff_artifact_dir / "diff_subplot",
        )

    def test_baseline_metadata_captures_dependency_versions(self):
        baseline_metadata_path = BASELINE_DIR / BASELINE_METADATA_FILENAME
        assert baseline_metadata_path.exists()

        metadata = json.loads(baseline_metadata_path.read_text(encoding="utf-8"))

        assert "python" in metadata
        assert "e3sm_diags_git_sha" in metadata
        assert "e3sm_unified_version" in metadata
        assert "e3sm_unified_feedstock_ref" in metadata
        assert "e3sm_unified_recipe_url" in metadata
        assert "e3sm_unified_recipe_version" in metadata

        for dependency in (
            "numpy",
            "pandas",
            "matplotlib",
            "cartopy",
            "xarray",
            "xcdat",
            "xesmf",
            "esmf",
            "esmpy",
            "xgcm",
        ):
            assert dependency in metadata
