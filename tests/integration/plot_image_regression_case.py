from __future__ import annotations

import os
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import numpy as np
import xarray as xr

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.plot.cosp_histogram_plot import plot as plot_cosp_histogram
from e3sm_diags.plot.lat_lon_plot import plot as plot_lat_lon
from e3sm_diags.plot.polar_plot import plot as plot_polar
from e3sm_diags.plot.zonal_mean_2d_plot import plot as plot_zonal_mean_2d

RenderFn = Callable[[str | Path], tuple[Path, ...]]

BASELINES_ROOT_DIR = Path("tests/integration/baselines")
SHARED_BASELINE_METADATA_PATH = BASELINES_ROOT_DIR / "baseline_metadata.json"


@dataclass(frozen=True)
class ImageRegressionCase:
    case_id: str
    baseline_dir: Path
    baseline_metadata_path: Path
    expected_image_filenames: tuple[str, ...]
    render: RenderFn
    mismatch_threshold: float = 0.0002


def _create_core_parameter(
    results_dir: str | Path,
    *,
    current_set: str,
    output_file: str,
    main_title: str,
) -> CoreParameter:
    parameter = CoreParameter()
    parameter.results_dir = str(results_dir)
    parameter.current_set = current_set
    parameter.case_id = "image_regression"
    parameter.output_file = output_file
    parameter.output_format = ["png"]
    parameter.output_format_subplot = ["png"]
    parameter.main_title = main_title
    parameter.test_name_yrs = "2000-2001"
    parameter.test_title = "Test"
    parameter.ref_name = "baseline"
    parameter.ref_name_yrs = "1990-1991"
    parameter.reference_title = "Reference"
    parameter.diff_name = "Test - Ref"
    parameter.diff_title = "Difference"
    parameter.test_colormap = (
        "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
    )
    parameter.reference_colormap = (
        "e3sm_diags/plot/colormaps/cet_rainbow_bgyrm_35_85_c71_r.rgb"
    )
    parameter.diff_colormap = "e3sm_diags/plot/colormaps/diverging_bwr.rgb"

    return parameter


def _get_output_dir(parameter: CoreParameter) -> Path:
    return Path(parameter.results_dir) / parameter.current_set / parameter.case_id


def _get_expected_output_paths(
    parameter: CoreParameter, expected_image_filenames: tuple[str, ...]
) -> tuple[Path, ...]:
    output_dir = _get_output_dir(parameter)
    return tuple(output_dir / filename for filename in expected_image_filenames)


def _create_latitudes(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(
        values,
        dims=("lat",),
        attrs={
            "axis": "Y",
            "long_name": "latitude",
            "standard_name": "latitude",
            "units": "degrees_north",
        },
    )


def _create_longitudes(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(
        values,
        dims=("lon",),
        attrs={
            "axis": "X",
            "long_name": "longitude",
            "standard_name": "longitude",
            "units": "degrees_east",
        },
    )


def _create_pressure_levels(values: np.ndarray) -> xr.DataArray:
    return xr.DataArray(
        values,
        dims=("plev",),
        attrs={
            "axis": "Z",
            "long_name": "pressure",
            "positive": "down",
            "standard_name": "air_pressure",
            "units": "mb",
        },
    )


def _create_2d_metrics_dict(
    da_test: xr.DataArray,
    da_ref: xr.DataArray,
    da_diff: xr.DataArray,
    *,
    unit_key: Literal["unit", "units"],
    units: str,
) -> MetricsDict:
    return {
        unit_key: units,
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


def _create_lat_lon_plot_data() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    lat = _create_latitudes(np.array([-60.0, -30.0, 0.0, 30.0, 60.0]))
    lon = _create_longitudes(
        np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
    )

    lon_grid, lat_grid = np.meshgrid(lon.data, lat.data)
    test_values = np.sin(np.deg2rad(lat_grid)) + 0.8 * np.cos(np.deg2rad(lon_grid))
    ref_values = 0.9 * np.sin(np.deg2rad(lat_grid + 5.0)) + 0.7 * np.cos(
        np.deg2rad(lon_grid - 15.0)
    )
    diff_values = test_values - ref_values

    da_test = xr.DataArray(
        test_values,
        name="test",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic field", "units": "K"},
    )
    da_ref = xr.DataArray(
        ref_values,
        name="ref",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic field", "units": "K"},
    )
    da_diff = xr.DataArray(
        diff_values,
        name="diff",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic field", "units": "K"},
    )

    return da_test, da_ref, da_diff


def _create_polar_plot_data() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    lat = _create_latitudes(np.array([55.0, 62.0, 69.0, 76.0, 83.0]))
    lon = _create_longitudes(
        np.array([0.0, 45.0, 90.0, 135.0, 180.0, 225.0, 270.0, 315.0])
    )

    lon_grid, lat_grid = np.meshgrid(lon.data, lat.data)
    test_values = (
        3.0 * np.sin(np.deg2rad(lat_grid - 50.0))
        + 1.5 * np.cos(np.deg2rad(lon_grid))
        + 0.4 * np.sin(np.deg2rad(lon_grid * 2.0))
    )
    ref_values = 2.7 * np.sin(np.deg2rad(lat_grid - 47.0)) + 1.2 * np.cos(
        np.deg2rad(lon_grid - 20.0)
    )
    diff_values = test_values - ref_values

    da_test = xr.DataArray(
        test_values,
        name="test",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic polar field", "units": "K"},
    )
    da_ref = xr.DataArray(
        ref_values,
        name="ref",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic polar field", "units": "K"},
    )
    da_diff = xr.DataArray(
        diff_values,
        name="diff",
        dims=("lat", "lon"),
        coords={"lat": lat.copy(), "lon": lon.copy()},
        attrs={"long_name": "synthetic polar field", "units": "K"},
    )

    return da_test, da_ref, da_diff


def _create_zonal_mean_2d_plot_data() -> tuple[
    xr.DataArray, xr.DataArray, xr.DataArray
]:
    lat = _create_latitudes(np.array([-90.0, -60.0, -30.0, 0.0, 30.0, 60.0, 90.0]))
    plev = _create_pressure_levels(
        np.array([1000.0, 850.0, 700.0, 500.0, 300.0, 200.0, 100.0])
    )

    lat_grid, plev_grid = np.meshgrid(lat.data, plev.data)
    test_values = (
        2.0 * np.cos(np.deg2rad(lat_grid))
        - 0.002 * plev_grid
        + 0.5 * np.sin(np.deg2rad(lat_grid * 3.0))
    )
    ref_values = (
        1.8 * np.cos(np.deg2rad(lat_grid + 8.0))
        - 0.0018 * plev_grid
        + 0.35 * np.sin(np.deg2rad(lat_grid * 2.0))
    )
    diff_values = test_values - ref_values

    da_test = xr.DataArray(
        test_values,
        name="test",
        dims=("plev", "lat"),
        coords={"plev": plev.copy(), "lat": lat.copy()},
        attrs={"long_name": "synthetic zonal mean field", "units": "K"},
    )
    da_ref = xr.DataArray(
        ref_values,
        name="ref",
        dims=("plev", "lat"),
        coords={"plev": plev.copy(), "lat": lat.copy()},
        attrs={"long_name": "synthetic zonal mean field", "units": "K"},
    )
    da_diff = xr.DataArray(
        diff_values,
        name="diff",
        dims=("plev", "lat"),
        coords={"plev": plev.copy(), "lat": lat.copy()},
        attrs={"long_name": "synthetic zonal mean field", "units": "K"},
    )

    return da_test, da_ref, da_diff


def _create_cosp_histogram_plot_data() -> tuple[
    xr.DataArray, xr.DataArray, xr.DataArray
]:
    tau = xr.DataArray(np.arange(6), dims=("tau",))
    ctp = xr.DataArray(np.arange(7), dims=("ctp",))

    test_values = np.array(
        [
            [0.5, 1.4, 2.1, 3.0, 4.4, 5.8],
            [0.9, 1.8, 3.2, 4.6, 6.1, 7.5],
            [1.6, 2.5, 4.1, 5.5, 7.2, 8.6],
            [2.2, 3.3, 4.9, 6.4, 8.0, 9.3],
            [3.0, 4.2, 5.8, 7.1, 8.9, 10.1],
            [3.6, 4.8, 6.5, 7.9, 9.6, 10.8],
            [4.1, 5.5, 7.0, 8.5, 10.1, 11.4],
        ]
    )
    ref_values = np.array(
        [
            [0.7, 1.2, 1.9, 2.8, 4.0, 5.2],
            [1.0, 1.7, 2.8, 4.1, 5.5, 6.8],
            [1.5, 2.3, 3.7, 5.0, 6.4, 7.9],
            [2.0, 3.0, 4.4, 5.8, 7.1, 8.4],
            [2.6, 3.7, 5.2, 6.6, 7.9, 9.1],
            [3.2, 4.3, 5.9, 7.2, 8.6, 9.8],
            [3.7, 4.9, 6.4, 7.8, 9.2, 10.4],
        ]
    )
    diff_values = test_values - ref_values

    da_test = xr.DataArray(
        test_values,
        name="CLDTOT_TEST",
        dims=("ctp", "tau"),
        coords={"ctp": ctp.copy(), "tau": tau.copy()},
        attrs={"long_name": "synthetic cosp histogram", "units": "%"},
    )
    da_ref = xr.DataArray(
        ref_values,
        name="CLDTOT_REF",
        dims=("ctp", "tau"),
        coords={"ctp": ctp.copy(), "tau": tau.copy()},
        attrs={"long_name": "synthetic cosp histogram", "units": "%"},
    )
    da_diff = xr.DataArray(
        diff_values,
        name="CLDTOT_DIFF",
        dims=("ctp", "tau"),
        coords={"ctp": ctp.copy(), "tau": tau.copy()},
        attrs={"long_name": "synthetic cosp histogram", "units": "%"},
    )

    return da_test, da_ref, da_diff


def render_lat_lon_plot_regression(results_dir: str | Path) -> tuple[Path, ...]:
    parameter = _create_core_parameter(
        results_dir,
        current_set="lat_lon",
        output_file="lat_lon_plot_regression",
        main_title="Synthetic Lat/Lon Regression",
    )
    parameter.contour_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    parameter.diff_levels = [-1.5, -0.75, 0.0, 0.75, 1.5]
    parameter.regions = ["45S45N-120E60W"]

    da_test, da_ref, da_diff = _create_lat_lon_plot_data()

    plot_lat_lon(
        parameter=parameter,
        da_test=da_test,
        da_ref=da_ref,
        da_diff=da_diff,
        metrics_dict=_create_2d_metrics_dict(
            da_test,
            da_ref,
            da_diff,
            unit_key="unit",
            units="K",
        ),
    )

    return _get_expected_output_paths(
        parameter,
        ("lat_lon_plot_regression.png", "lat_lon_plot_regression.2.png"),
    )


def render_polar_plot_regression(results_dir: str | Path) -> tuple[Path, ...]:
    parameter = _create_core_parameter(
        results_dir,
        current_set="polar",
        output_file="polar_plot_regression",
        main_title="Synthetic Polar Regression",
    )
    parameter.contour_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]
    parameter.diff_levels = [-1.5, -0.75, 0.0, 0.75, 1.5]
    parameter.var_region = "60N90N"

    da_test, da_ref, da_diff = _create_polar_plot_data()

    plot_polar(
        parameter=parameter,
        da_test=da_test,
        da_ref=da_ref,
        da_diff=da_diff,
        metrics_dict=_create_2d_metrics_dict(
            da_test,
            da_ref,
            da_diff,
            unit_key="units",
            units="K",
        ),
    )

    return _get_expected_output_paths(
        parameter,
        ("polar_plot_regression.png", "polar_plot_regression.2.png"),
    )


def render_zonal_mean_2d_plot_regression(results_dir: str | Path) -> tuple[Path, ...]:
    parameter = _create_core_parameter(
        results_dir,
        current_set="zonal_mean_2d",
        output_file="zonal_mean_2d_plot_regression",
        main_title="Synthetic Zonal Mean 2D Regression",
    )
    parameter.contour_levels = [-1.0, 0.0, 1.0, 2.0, 3.0]
    parameter.diff_levels = [-0.8, -0.4, 0.0, 0.4, 0.8]
    parameter.plot_plevs = True
    parameter.plevs = [1000, 850, 700, 500, 300, 200, 100]

    da_test, da_ref, da_diff = _create_zonal_mean_2d_plot_data()

    plot_zonal_mean_2d(
        parameter=parameter,
        da_test=da_test,
        da_ref=da_ref,
        da_diff=da_diff,
        metrics_dict=_create_2d_metrics_dict(
            da_test,
            da_ref,
            da_diff,
            unit_key="units",
            units="K",
        ),
    )

    return _get_expected_output_paths(
        parameter,
        ("zonal_mean_2d_plot_regression.png", "zonal_mean_2d_plot_regression.2.png"),
    )


def render_cosp_histogram_plot_regression(results_dir: str | Path) -> tuple[Path, ...]:
    parameter = _create_core_parameter(
        results_dir,
        current_set="cosp_histogram",
        output_file="cosp_histogram_plot_regression",
        main_title="Synthetic COSP Histogram Regression",
    )
    parameter.contour_levels = [1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0]
    parameter.diff_levels = [-2.0, -1.0, 0.0, 1.0, 2.0]

    da_test, da_ref, da_diff = _create_cosp_histogram_plot_data()

    plot_cosp_histogram(
        parameter=parameter,
        da_test=da_test,
        da_ref=da_ref,
        da_diff=da_diff,
    )

    return _get_expected_output_paths(
        parameter,
        ("cosp_histogram_plot_regression.png", "cosp_histogram_plot_regression.2.png"),
    )


IMAGE_REGRESSION_CASES = (
    ImageRegressionCase(
        case_id="lat_lon",
        baseline_dir=BASELINES_ROOT_DIR / "lat_lon_plot",
        baseline_metadata_path=SHARED_BASELINE_METADATA_PATH,
        expected_image_filenames=(
            "lat_lon_plot_regression.png",
            "lat_lon_plot_regression.2.png",
        ),
        render=render_lat_lon_plot_regression,
    ),
    ImageRegressionCase(
        case_id="polar",
        baseline_dir=BASELINES_ROOT_DIR / "polar_plot",
        baseline_metadata_path=SHARED_BASELINE_METADATA_PATH,
        expected_image_filenames=(
            "polar_plot_regression.png",
            "polar_plot_regression.2.png",
        ),
        render=render_polar_plot_regression,
        mismatch_threshold=0.005,
    ),
    ImageRegressionCase(
        case_id="zonal_mean_2d",
        baseline_dir=BASELINES_ROOT_DIR / "zonal_mean_2d_plot",
        baseline_metadata_path=SHARED_BASELINE_METADATA_PATH,
        expected_image_filenames=(
            "zonal_mean_2d_plot_regression.png",
            "zonal_mean_2d_plot_regression.2.png",
        ),
        render=render_zonal_mean_2d_plot_regression,
    ),
    ImageRegressionCase(
        case_id="cosp_histogram",
        baseline_dir=BASELINES_ROOT_DIR / "cosp_histogram_plot",
        baseline_metadata_path=SHARED_BASELINE_METADATA_PATH,
        expected_image_filenames=(
            "cosp_histogram_plot_regression.png",
            "cosp_histogram_plot_regression.2.png",
        ),
        render=render_cosp_histogram_plot_regression,
    ),
)
IMAGE_REGRESSION_CASES_BY_ID = {case.case_id: case for case in IMAGE_REGRESSION_CASES}


def get_image_regression_artifact_dir(tmp_path: Path) -> Path:
    base_dir = os.environ.get("IMAGE_REGRESSION_ARTIFACT_DIR")

    if base_dir is None:
        return tmp_path / "image_regression_artifacts"

    return Path(base_dir)
