from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.integration.image_regression import assert_image_matches_baseline
from tests.integration.plot_image_regression_case import (
    IMAGE_REGRESSION_CASES,
    ImageRegressionCase,
    get_image_regression_artifact_dir,
)


@pytest.mark.image_regression
class TestPlotImageRegressions:
    @pytest.mark.parametrize(
        "case",
        IMAGE_REGRESSION_CASES,
        ids=[case.case_id for case in IMAGE_REGRESSION_CASES],
    )
    def test_plot_matches_expected_baselines(
        self,
        tmp_path: Path,
        case: ImageRegressionCase,
    ):
        generated_images = case.render(tmp_path / "results")
        diff_artifact_dir = get_image_regression_artifact_dir(tmp_path) / case.case_id

        for generated_image, image_filename in zip(
            generated_images, case.expected_image_filenames, strict=True
        ):
            runtime_metadata_path = assert_image_matches_baseline(
                image_name=image_filename,
                path_to_actual_png=generated_image,
                path_to_expected_png=case.baseline_dir / image_filename,
                artifact_dir=diff_artifact_dir / Path(image_filename).stem,
                mismatch_threshold=case.get_mismatch_threshold(),
            )
            assert runtime_metadata_path.exists()

    @pytest.mark.parametrize(
        "case",
        IMAGE_REGRESSION_CASES,
        ids=[f"{case.case_id}_metadata" for case in IMAGE_REGRESSION_CASES],
    )
    def test_baseline_metadata_captures_dependency_versions(
        self,
        case: ImageRegressionCase,
    ):
        baseline_metadata_path = case.baseline_metadata_path
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
