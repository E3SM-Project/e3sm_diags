import numpy as np
import pytest
import xarray as xr
from xarray.testing import assert_allclose

from e3sm_diags.metrics.metrics import correlation, get_weights, rmse, spatial_avg, std


class TestGetWeights:
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

    def test_returns_weights_for_x_y_axes(self):
        expected = xr.DataArray(
            name="lat_lon_wts",
            data=np.array(
                [[0.01745241, 0.01744709], [0.01745241, 0.01744709]], dtype="float64"
            ),
            coords={"lon": self.ds.lon, "lat": self.ds.lat},
        )
        result = get_weights(self.ds)

        assert_allclose(expected, result)


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

    def test_returns_spatial_avg_for_x_y(self):
        expected = [1.5, 1.333299, 1.5]
        result = spatial_avg(self.ds, "ts")

        np.testing.assert_allclose(expected, result, atol=1e-5, rtol=1e-5)

    def test_returns_spatial_avg_for_x_y_as_xr_dataarray(self):
        expected = [1.5, 1.333299, 1.5]
        result = spatial_avg(self.ds, "ts", as_list=False)

        assert isinstance(result, xr.DataArray)
        np.testing.assert_allclose(expected, result, atol=1e-5, rtol=1e-5)


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

    def test_returns_weighted_std_for_x_y_axes(self):
        expected = [0.5, 0.47139255, 0.5]
        result = std(self.ds, "ts")

        np.testing.assert_allclose(expected, result)


class TestCorrelation:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.var_key = "ts"
        self.ds_a = xr.Dataset(
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

        self.ds_a[self.var_key] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds_a.lat, "lon": self.ds_a.lon, "time": self.ds_a.time},
            dims=["time", "lat", "lon"],
        )

        # Bounds are used to generate weights.
        self.ds_a["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds_a["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

        self.ds_b = self.ds_a.copy()
        self.ds_b[self.var_key] = xr.DataArray(
            data=np.array(
                [
                    [[1, 2.25], [0.925, 2.10]],
                    [[np.nan, 1.2], [1.1, 2]],
                    [[2, 1.1], [1.1, 2]],
                ]
            ),
            coords={"lat": self.ds_a.lat, "lon": self.ds_a.lon, "time": self.ds_a.time},
            dims=["time", "lat", "lon"],
        )

    def test_returns_weighted_correlation_on_x_y_axes(self):
        expected = [0.99525143, 0.99484914, 1]

        result = correlation(self.ds_a, self.ds_b, self.var_key)

        np.testing.assert_allclose(expected, result)


class TestRmse:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.var_key = "ts"
        self.ds_a = xr.Dataset(
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

        self.ds_a[self.var_key] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds_a.lat, "lon": self.ds_a.lon, "time": self.ds_a.time},
            dims=["time", "lat", "lon"],
        )

        # Bounds are used to generate weights.
        self.ds_a["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds_a["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

        self.ds_b = self.ds_a.copy()
        self.ds_b[self.var_key] = xr.DataArray(
            data=np.array(
                [
                    [[1, 2.25], [0.925, 2.10]],
                    [[np.nan, 1.2], [1.1, 2]],
                    [[2, 1.1], [1.1, 2]],
                ]
            ),
            coords={"lat": self.ds_a.lat, "lon": self.ds_a.lon, "time": self.ds_a.time},
            dims=["time", "lat", "lon"],
        )

    def test_returns_weighted_rmse_on_x_y_axes(self):
        expected = [0.13976063, 0.12910862, 0.07071068]

        result = rmse(self.ds_a, self.ds_b, self.var_key)

        np.testing.assert_allclose(expected, result)
