import numpy as np
import pytest
import xarray as xr
from xarray.testing import assert_allclose

from e3sm_diags.metrics.metrics import (  # noqa: F401
    correlation,
    get_weights,
    rmse,
    spatial_avg,
    std,
)


class TestGetWeights:
    pass


class TestSpatialAvg:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.ds = xr.Dataset(
            coords={
                "lat": xr.DataArray(
                    data=[0, 1], dims="lat", attrs={"bounds": "lat_bnds", "axis": "Y"}
                ),
                "lon": xr.DataArray(
                    data=[0, 1], dims="lon", attrs={"bounds": "lon_bnds", "axis": "X"}
                ),
                "time": xr.DataArray(data=[1, 2, 3], dims="time"),
            },
        )

        self.ds["ts"] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )

        # Bounds are used to generate weights.
        self.ds["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

    def test_raises_error_with_invalid_axis_arg(self):
        with pytest.raises(ValueError):
            spatial_avg(self.ds, "ts", axis=["T"])

        with pytest.raises(ValueError):
            spatial_avg(self.ds, "ts", axis="T")  # type: ignore

    def test_returns_spatial_avg_for_x_y_axes(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([1.5, 1.3333, 1.5]),
            coords={"time": self.ds.time},
            dims=["time"],
        )
        result = spatial_avg(self.ds, "ts")

        assert_allclose(expected, result)

    def test_returns_spatial_avg_for_x_axis(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.5, 1.5], [1, 1.5], [1.5, 1.5]]),
            coords={"time": self.ds.time, "lat": self.ds.lat},
            dims=["time", "lat"],
        )
        result = spatial_avg(self.ds, "ts", axis=["X"])

        assert_allclose(expected, result)

    def test_returns_spatial_avg_for_y_axis(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0, 2.0], [1, 1.499924], [1.50008, 1.499924]]),
            coords={"time": self.ds.time, "lon": self.ds.lon},
            dims=["time", "lon"],
        )
        result = spatial_avg(self.ds, "ts", axis=["Y"])

        assert_allclose(expected, result)


class TestStd:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.ds = xr.Dataset(
            coords={
                "lat": xr.DataArray(
                    data=[0, 1], dims="lat", attrs={"bounds": "lat_bnds", "axis": "Y"}
                ),
                "lon": xr.DataArray(
                    data=[0, 1], dims="lon", attrs={"bounds": "lon_bnds", "axis": "X"}
                ),
                "time": xr.DataArray(data=[1, 2, 3], dims="time"),
            },
        )

        self.ds["ts"] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )

        # Bounds are used to generate weights.
        self.ds["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

    def test_raises_error_with_invalid_axis_arg(self):
        with pytest.raises(ValueError):
            std(self.ds, "ts", axis=["T"])

        with pytest.raises(ValueError):
            std(self.ds, "ts", axis="T")

    def test_returns_weighted_std_for_x_y_axes(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([0.5, 0.47139255, 0.5]),
            coords={"time": self.ds.time},
            dims=["time"],
        )
        result = std(self.ds, "ts")

        assert_allclose(expected, result)

    def test_returns_weighted_std_for_x_axis(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([[0.5, 0.5], [0.0, 0.5], [0.5, 0.5]]),
            coords={"time": self.ds.time, "lat": self.ds.lat},
            dims=["time", "lat"],
        )
        result = std(self.ds, "ts", axis=["X"])

        assert_allclose(expected, result)

    def test_returns_weighted_std_for_y_axis(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([[0.0, 0.0], [0.0, 0.5], [0.5, 0.5]]),
            coords={"time": self.ds.time, "lon": self.ds.lon},
            dims=["time", "lon"],
        )
        result = std(self.ds, "ts", axis=["Y"])

        assert_allclose(expected, result)


class TestCorrelation:
    pass


class TestRMSE:
    pass
