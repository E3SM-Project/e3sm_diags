from __future__ import annotations

import os
from collections.abc import Callable, Mapping
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
COMPAT_TOLERANCE_PROFILE = "e3sm_unified_latest_released_compat"


@dataclass(frozen=True)
class ImageRegressionCase:
    case_id: str
    baseline_dir: Path
    baseline_metadata_path: Path
    expected_image_filenames: tuple[str, ...]
    render: RenderFn
    mismatch_threshold: float = 0.0002
    compat_mismatch_threshold: float | None = None
    compat_mismatch_thresholds_by_image: Mapping[str, float] | None = None

    def get_mismatch_threshold(self, image_filename: str) -> float:
        if (
            os.environ.get("IMAGE_REGRESSION_TOLERANCE_PROFILE")
            == COMPAT_TOLERANCE_PROFILE
        ):
            if self.compat_mismatch_thresholds_by_image is not None:
                threshold = self.compat_mismatch_thresholds_by_image.get(image_filename)
                if threshold is not None:
                    return threshold

            if self.compat_mismatch_threshold is not None:
                return self.compat_mismatch_threshold

        return self.mismatch_threshold


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


def _wrapped_lon_delta(
    lon_values: np.ndarray, lon_center: float
) -> np.ndarray:
    """Get wrapped longitude deltas in degrees."""
    return ((lon_values - lon_center + 180.0) % 360.0) - 180.0


def _gaussian_lat_lon(
    lat_grid: np.ndarray,
    lon_grid: np.ndarray,
    *,
    lat_center: float,
    lon_center: float,
    lat_sigma: float,
    lon_sigma: float,
) -> np.ndarray:
    """Get deterministic 2D Gaussian on lat/lon grid."""
    lon_delta = _wrapped_lon_delta(lon_grid, lon_center)
    return np.exp(
        -0.5
        * (
            ((lat_grid - lat_center) / lat_sigma) ** 2
            + (lon_delta / lon_sigma) ** 2
        )
    )


def _gaussian_lat_plev(
    lat_grid: np.ndarray,
    plev_grid: np.ndarray,
    *,
    lat_center: float,
    lat_sigma: float,
    plev_center: float,
    log_plev_sigma: float,
) -> np.ndarray:
    """Get deterministic 2D Gaussian on latitude-pressure grid."""
    log_plev_delta = np.log(plev_grid / plev_center)
    return np.exp(
        -0.5
        * (
            ((lat_grid - lat_center) / lat_sigma) ** 2
            + (log_plev_delta / log_plev_sigma) ** 2
        )
    )


def _create_lat_lon_plot_data() -> tuple[xr.DataArray, xr.DataArray, xr.DataArray]:
    lat = _create_latitudes(np.arange(-60.0, 65.0, 5.0))
    lon = _create_longitudes(np.arange(0.0, 360.0, 5.0))

    lon_grid, lat_grid = np.meshgrid(lon.data, lat.data)

    background_test = (
        0.22 * np.cos(np.deg2rad(lat_grid * 1.2))
        + 0.20
        * np.sin(np.deg2rad(lon_grid - 118.0))
        * np.cos(np.deg2rad(lat_grid * 0.8))
        - 0.14 * (lat_grid / 20.0)
    )
    background_ref = (
        0.20 * np.cos(np.deg2rad((lat_grid + 2.0) * 1.15))
        + 0.18
        * np.sin(np.deg2rad(lon_grid - 125.0))
        * np.cos(np.deg2rad((lat_grid + 1.0) * 0.8))
        - 0.12 * ((lat_grid + 1.0) / 20.0)
    )

    test_values = (
        background_test
        + 1.02
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=8.0,
            lon_center=145.0,
            lat_sigma=11.0,
            lon_sigma=18.0,
        )
        - 0.46
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=-12.0,
            lon_center=168.0,
            lat_sigma=9.0,
            lon_sigma=15.0,
        )
        + 0.28
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=20.0,
            lon_center=112.0,
            lat_sigma=10.0,
            lon_sigma=13.0,
        )
        - 0.56
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=-8.0,
            lon_center=176.0,
            lat_sigma=8.0,
            lon_sigma=11.0,
        )
    )
    ref_values = (
        background_ref
        + 0.90
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=10.0,
            lon_center=150.0,
            lat_sigma=12.0,
            lon_sigma=20.0,
        )
        - 0.34
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=-9.0,
            lon_center=162.0,
            lat_sigma=10.0,
            lon_sigma=17.0,
        )
        + 0.20
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=18.0,
            lon_center=120.0,
            lat_sigma=11.0,
            lon_sigma=16.0,
        )
        - 0.42
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=-6.0,
            lon_center=172.0,
            lat_sigma=9.0,
            lon_sigma=13.0,
        )
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
    lat = _create_latitudes(np.arange(50.0, 91.0, 4.0))
    lon = _create_longitudes(np.arange(0.0, 360.0, 15.0))

    lon_grid, lat_grid = np.meshgrid(lon.data, lat.data)
    radial_distance = (90.0 - lat_grid) / 40.0

    test_values = (
        1.32
        - 1.72 * radial_distance
        + 0.42
        * np.cos(np.deg2rad(lon_grid - 40.0))
        * (1.0 - 0.35 * radial_distance)
        + 0.18 * np.sin(np.deg2rad(2.0 * lon_grid + 20.0))
        + 0.74
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=78.0,
            lon_center=110.0,
            lat_sigma=6.0,
            lon_sigma=30.0,
        )
        - 0.58
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=67.0,
            lon_center=248.0,
            lat_sigma=7.0,
            lon_sigma=38.0,
        )
    )
    ref_values = (
        1.22
        - 1.60 * radial_distance
        + 0.35
        * np.cos(np.deg2rad(lon_grid - 58.0))
        * (1.0 - 0.30 * radial_distance)
        + 0.12 * np.sin(np.deg2rad(2.0 * lon_grid - 8.0))
        + 0.61
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=76.0,
            lon_center=128.0,
            lat_sigma=7.0,
            lon_sigma=34.0,
        )
        - 0.50
        * _gaussian_lat_lon(
            lat_grid,
            lon_grid,
            lat_center=69.0,
            lon_center=232.0,
            lat_sigma=8.0,
            lon_sigma=42.0,
        )
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
    lat = _create_latitudes(np.linspace(-90.0, 90.0, 25))
    plev = _create_pressure_levels(
        np.array(
            [
                1000.0,
                925.0,
                850.0,
                775.0,
                700.0,
                600.0,
                500.0,
                400.0,
                300.0,
                250.0,
                200.0,
                150.0,
                100.0,
                70.0,
            ]
        )
    )

    lat_grid, plev_grid = np.meshgrid(lat.data, plev.data)

    test_values = (
        0.34 * np.cos(np.deg2rad(lat_grid))
        - 0.92 * (plev_grid / 1000.0) ** 0.85 * np.sin(np.deg2rad(np.abs(lat_grid))) ** 1.1
        + 2.20
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=5.0,
            lat_sigma=24.0,
            plev_center=220.0,
            log_plev_sigma=0.48,
        )
        + 0.58
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=42.0,
            lat_sigma=17.0,
            plev_center=520.0,
            log_plev_sigma=0.40,
        )
        + 0.42
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=-46.0,
            lat_sigma=19.0,
            plev_center=580.0,
            log_plev_sigma=0.43,
        )
        + 0.10
        * np.sin(np.deg2rad(lat_grid * 2.0))
        * np.exp(-0.5 * (np.log(plev_grid / 650.0) / 0.55) ** 2)
    )
    ref_values = (
        0.30 * np.cos(np.deg2rad(lat_grid + 4.0))
        - 0.86 * (plev_grid / 1000.0) ** 0.82 * np.sin(np.deg2rad(np.abs(lat_grid + 4.0))) ** 1.05
        + 2.00
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=8.0,
            lat_sigma=26.0,
            plev_center=240.0,
            log_plev_sigma=0.50,
        )
        + 0.50
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=40.0,
            lat_sigma=18.0,
            plev_center=540.0,
            log_plev_sigma=0.42,
        )
        + 0.48
        * _gaussian_lat_plev(
            lat_grid,
            plev_grid,
            lat_center=-42.0,
            lat_sigma=20.0,
            plev_center=610.0,
            log_plev_sigma=0.45,
        )
        + 0.08
        * np.sin(np.deg2rad((lat_grid + 6.0) * 2.0))
        * np.exp(-0.5 * (np.log(plev_grid / 620.0) / 0.58) ** 2)
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
            [1.8, 2.3, 3.0, 2.6, 1.9, 1.5],
            [2.1, 2.8, 3.8, 3.4, 2.4, 1.8],
            [1.2, 1.7, 4.2, 5.1, 3.6, 2.9],
            [0.9, 1.4, 4.8, 5.8, 4.1, 3.2],
            [0.7, 1.1, 2.6, 3.3, 4.5, 5.2],
            [0.5, 0.9, 2.0, 2.7, 5.1, 6.0],
            [0.3, 0.6, 1.4, 2.0, 4.8, 5.6],
        ]
    )
    diff_values = np.array(
        [
            [0.4, 0.3, -0.2, -0.3, -0.4, -0.2],
            [0.3, 0.2, -0.1, -0.2, -0.3, -0.2],
            [0.2, 0.1, 0.4, 0.5, -0.2, -0.4],
            [0.1, 0.0, 0.5, 0.6, -0.3, -0.5],
            [-0.1, -0.1, 0.2, 0.3, 0.5, 0.6],
            [-0.1, -0.2, 0.1, 0.2, 0.7, 0.8],
            [-0.1, -0.1, 0.0, 0.1, 0.6, 0.7],
        ]
    )
    ref_values = test_values - diff_values

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
    parameter.regions = ["W_Pacific"]

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
        # The latest released E3SM-Unified stack shows renderer-only drift in
        # polar linework and clipping. Linux compat runs have measured roughly
        # 3.6% mismatched pixels for the full figure and 8.2% for the cropped
        # difference subplot while the filled field remains visually
        # equivalent, so keep the strict default threshold and relax only the
        # compat profile.
        compat_mismatch_threshold=0.04,
        compat_mismatch_thresholds_by_image={
            "polar_plot_regression.2.png": 0.09,
        },
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
