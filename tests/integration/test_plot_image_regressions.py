from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.lat_lon_plot import plot
from tests.integration.image_regression import (
    BASELINE_METADATA_FILENAME,
    assert_image_matches_baseline,
)

BASELINE_DIR = Path("tests/integration/baselines/lat_lon_plot")


@pytest.mark.image_regression
class TestLatLonPlotImageRegressions:
    def test_lat_lon_plot_matches_expected_baselines(self, tmp_path: Path):
        parameter = self._create_parameter(tmp_path)
        da_test, da_ref, da_diff = self._create_plot_data()

        plot(
            parameter=parameter,
            da_test=da_test,
            da_ref=da_ref,
            da_diff=da_diff,
            metrics_dict=self._create_metrics_dict(da_test, da_ref, da_diff),
        )

        output_dir = Path(parameter.results_dir) / parameter.current_set / parameter.case_id
        diff_artifact_dir = tmp_path / "image_regression_artifacts"

        runtime_metadata_path = assert_image_matches_baseline(
            image_name="lat_lon_plot_regression.png",
            path_to_actual_png=output_dir / f"{parameter.output_file}.png",
            path_to_expected_png=BASELINE_DIR / "lat_lon_plot_regression.png",
            artifact_dir=diff_artifact_dir / "main_plot",
        )
        assert runtime_metadata_path.exists()

        assert_image_matches_baseline(
            image_name="lat_lon_plot_regression.2.png",
            path_to_actual_png=output_dir / f"{parameter.output_file}.2.png",
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

    def _create_parameter(self, tmp_path: Path) -> CoreParameter:
        parameter = CoreParameter()
        parameter.results_dir = str(tmp_path / "results")
        parameter.current_set = "lat_lon"
        parameter.case_id = "image_regression"
        parameter.output_file = "lat_lon_plot_regression"
        parameter.output_format = ["png"]
        parameter.output_format_subplot = ["png"]
        parameter.main_title = "Synthetic Lat/Lon Regression"
        parameter.test_name_yrs = "2000-2001"
        parameter.test_title = "Test"
        parameter.ref_name = "baseline"
        parameter.ref_name_yrs = "1990-1991"
        parameter.reference_title = "Reference"
        parameter.diff_title = "Difference"
        parameter.test_colormap = "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
        parameter.reference_colormap = (
            "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
        )
        parameter.diff_colormap = "e3sm_diags/plot/colormaps/diverging_bwr.rgb"
        parameter.contour_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
        parameter.diff_levels = [-1.5, -0.75, 0.0, 0.75, 1.5]
        parameter.regions = ["45S45N-120E60W"]

        return parameter

    def _create_plot_data(self) -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
        lat = xr.DataArray(
            np.array([-60.0, -30.0, 0.0, 30.0, 60.0]),
            dims=("lat",),
            attrs={
                "axis": "Y",
                "long_name": "latitude",
                "standard_name": "latitude",
                "units": "degrees_north",
            },
        )
        lon = xr.DataArray(
            np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0]),
            dims=("lon",),
            attrs={
                "axis": "X",
                "long_name": "longitude",
                "standard_name": "longitude",
                "units": "degrees_east",
            },
        )

        lon_grid, lat_grid = np.meshgrid(lon.data, lat.data)
        base = np.sin(np.deg2rad(lat_grid)) + 0.8 * np.cos(np.deg2rad(lon_grid))
        ref = 0.9 * np.sin(np.deg2rad(lat_grid + 5.0)) + 0.7 * np.cos(
            np.deg2rad(lon_grid - 15.0)
        )
        diff = base - ref

        data_arrays = []
        for name, values in (("test", base), ("ref", ref), ("diff", diff)):
            data_arrays.append(
                xr.DataArray(
                    values,
                    name=name,
                    dims=("lat", "lon"),
                    coords={"lat": lat.copy(), "lon": lon.copy()},
                    attrs={"long_name": "synthetic field", "units": "K"},
                )
            )

        return tuple(data_arrays)  # type: ignore[return-value]

    def _create_metrics_dict(
        self, da_test: xr.DataArray, da_ref: xr.DataArray, da_diff: xr.DataArray
    ) -> dict[str, object]:
        return {
            "unit": "K",
            "test": {
                "min": float(da_test.min()),
                "mean": float(da_test.mean()),
                "max": float(da_test.max()),
            },
            "ref": {
                "min": float(da_ref.min()),
                "mean": float(da_ref.mean()),
                "max": float(da_ref.max()),
            },
            "diff": {
                "min": float(da_diff.min()),
                "mean": float(da_diff.mean()),
                "max": float(da_diff.max()),
            },
            "misc": {
                "rmse": float(np.sqrt(np.mean(np.square(da_diff.data)))),
                "corr": float(np.corrcoef(da_test.data.ravel(), da_ref.data.ravel())[0, 1]),
            },
        }
