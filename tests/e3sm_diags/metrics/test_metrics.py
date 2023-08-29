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

    def test_returns_weights_for_x_axis(self):
        expected = xr.DataArray(
            name="lon_wts",
            data=np.array([1, 1], dtype="float64"),
            coords={"lon": self.ds.lon},
        )
        result = get_weights(self.ds)

        assert_allclose(expected, result)

    def test_returns_weights_for_y_axis(self):
        expected = xr.DataArray(
            name="lat_wts",
            data=np.array([0.01745241, 0.01744709], dtype="float64"),
            coords={"lat": self.ds.lat},
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

    def test_returns_spatial_avg_for_x_y_axes(self):
        expected = xr.DataArray(
            name="ts",
            data=np.array([1.5, 1.3333, 1.5]),
            coords={"time": self.ds.time},
            dims=["time"],
        )
        result = spatial_avg(self.ds, "ts")

        assert_allclose(expected, result)

    def test_returns_serialized_spatial_avg_for_x_y(self):
        expected = [1.5, 1.3333, 1.5]
        result = spatial_avg(self.ds, "ts", serialize=True)

        np.testing.assert_allclose(expected, result)


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
        expected = xr.DataArray(
            name="ts",
            data=np.array([0.5, 0.47139255, 0.5]),
            coords={"time": self.ds.time},
            dims=["time"],
        )
        result = std(self.ds, "ts")

        assert_allclose(expected, result)

    def test_returns_serialized_weighted_std_for_x_y_axes(self):
        expected = [0.5, 0.47139255, 0.5]
        result = std(self.ds, "ts", serialize=True)

        np.testing.assert_allclose(expected, result)


class TestCorrelation:
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

        self.ds["ts_model"] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )
        self.ds["ts_obs"] = xr.DataArray(
            data=np.array(
                [
                    [[1, 2.25], [0.925, 2.10]],
                    [[np.nan, 1.2], [1.1, 2]],
                    [[2, 1.1], [1.1, 2]],
                ]
            ),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )
        # Bounds are used to generate weights.
        self.ds["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

    def test_returns_weighted_correlation_on_x_y_axes(self):
        expected = xr.DataArray(
            data=np.array([0.99525143, np.nan, 1], dtype="float64"),
            coords={"time": self.ds.time},
        )

        weights = get_weights(self.ds)
        result = correlation(self.ds.ts_model, self.ds.ts_obs, weights=weights)

        assert_allclose(expected, result)

    def test_returns_serialized_weighted_correlation_on_x_y_axes(self):
        expected = [0.99525143, np.nan, 1]

        weights = get_weights(self.ds)
        result = correlation(self.ds.ts_model, self.ds.ts_obs, weights=weights)

        np.testing.assert_allclose(expected, result)


class TestRmse:
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

        self.ds["ts_model"] = xr.DataArray(
            data=np.array([[[1, 2], [1, 2]], [[np.nan, 1], [1, 2]], [[2, 1], [1, 2]]]),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )
        self.ds["ts_obs"] = xr.DataArray(
            data=np.array(
                [
                    [[1, 2.25], [0.925, 2.10]],
                    [[np.nan, 1.2], [1.1, 2]],
                    [[2, 1.1], [1.1, 2]],
                ]
            ),
            coords={"lat": self.ds.lat, "lon": self.ds.lon, "time": self.ds.time},
            dims=["time", "lat", "lon"],
        )
        # Bounds are used to generate weights.
        self.ds["lat_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lat", "bnds"])
        self.ds["lon_bnds"] = xr.DataArray([[0, 1], [1, 2]], dims=["lon", "bnds"])

    def test_returns_weighted_rmse_on_x_y_axes(self):
        expected = xr.DataArray(
            data=np.array([0.13976063, np.nan, 0.07071068], dtype="float64"),
            coords={"time": self.ds.time},
        )

        weights = get_weights(self.ds)
        result = rmse(self.ds.ts_model, self.ds.ts_obs, weights=weights)

        assert_allclose(expected, result)

    def test_returns_serialized_weighted_rmse_on_x_y_axes(self):
        expected = [0.13976063, np.nan, 0.07071068]

        weights = get_weights(self.ds)
        result = rmse(
            self.ds.ts_model,
            self.ds.ts_obs,
            weights=weights,
            serialize=True,
        )

        np.testing.assert_allclose(expected, result)
