from __future__ import annotations

import os
from pathlib import Path

import numpy as np
import xarray as xr

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.lat_lon_plot import plot

BASELINE_DIR = Path("tests/integration/baselines/lat_lon_plot")
BASELINE_IMAGE_FILENAMES = (
    "lat_lon_plot_regression.png",
    "lat_lon_plot_regression.2.png",
)


def create_lat_lon_parameter(results_dir: str | Path) -> CoreParameter:
    parameter = CoreParameter()
    parameter.results_dir = str(results_dir)
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
    parameter.test_colormap = (
        "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
    )
    parameter.reference_colormap = (
        "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
    )
    parameter.diff_colormap = "e3sm_diags/plot/colormaps/diverging_bwr.rgb"
    parameter.contour_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    parameter.diff_levels = [-1.5, -0.75, 0.0, 0.75, 1.5]
    parameter.regions = ["45S45N-120E60W"]

    return parameter


def create_lat_lon_plot_data() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
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

    data_arrays: list[xr.DataArray] = []
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

    return tuple(data_arrays)


def create_metrics_dict(
    da_test: xr.DataArray, da_ref: xr.DataArray, da_diff: xr.DataArray
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


def render_lat_lon_plot_regression(
    results_dir: str | Path,
) -> tuple[CoreParameter, tuple[Path, ...]]:
    parameter = create_lat_lon_parameter(results_dir)
    da_test, da_ref, da_diff = create_lat_lon_plot_data()

    plot(
        parameter=parameter,
        da_test=da_test,
        da_ref=da_ref,
        da_diff=da_diff,
        metrics_dict=create_metrics_dict(da_test, da_ref, da_diff),
    )

    output_dir = Path(parameter.results_dir) / parameter.current_set / parameter.case_id
    generated_images = tuple(
        output_dir / filename for filename in BASELINE_IMAGE_FILENAMES
    )

    return parameter, generated_images


def get_image_regression_artifact_dir(tmp_path: Path) -> Path:
    base_dir = os.environ.get("IMAGE_REGRESSION_ARTIFACT_DIR")

    if base_dir is None:
        return tmp_path / "image_regression_artifacts"

    return Path(base_dir)
