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

    def test_raises_error_with_invalid_axis_arg(self):
        with pytest.raises(ValueError):
            get_weights(self.ds, axis=["T"])

        with pytest.raises(ValueError):
            get_weights(self.ds, axis="T")  # type: ignore

    def test_returns_weights_for_x_y_axes(self):
        expected = xr.DataArray(
            name="lat_lon_wts",
            data=np.array(
                [[0.01745241, 0.01744709], [0.01745241, 0.01744709]], dtype="float64"
            ),
            coords={"lon": self.ds.lon, "lat": self.ds.lat},
        )
        result = get_weights(self.ds, axis=["X", "Y"])

        assert_allclose(expected, result)

    def test_returns_weights_for_x_axis(self):
        expected = xr.DataArray(
            name="lon_wts",
            data=np.array([1, 1], dtype="float64"),
            coords={"lon": self.ds.lon},
        )
        result = get_weights(self.ds, axis=["X"])

        assert_allclose(expected, result)

    def test_returns_weights_for_y_axis(self):
        expected = xr.DataArray(
            name="lat_wts",
            data=np.array([0.01745241, 0.01744709], dtype="float64"),
            coords={"lat": self.ds.lat},
        )
        result = get_weights(self.ds, axis=["Y"])

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

        weights = get_weights(self.ds, axis=["X", "Y"])
        result = correlation(
            self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["X", "Y"]
        )

        assert_allclose(expected, result)

    def test_returns_weighted_correlation_on_x_axis(self):
        expected = xr.DataArray(
            data=np.array([[1.0, 1.0], [np.nan, 1.0], [1.0, 1.0]], dtype="float64"),
            coords={"time": self.ds.time, "lat": self.ds.lat},
        )

        weights = get_weights(self.ds, axis=["X"])
        result = correlation(
            self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["X"]
        )

        assert_allclose(expected, result)

    def test_returns_weighted_correlation_on_y_axis(self):
        expected = xr.DataArray(
            data=np.array(
                [[np.nan, np.nan], [np.nan, 1.0], [1.0, 1.0]], dtype="float64"
            ),
            coords={"time": self.ds.time, "lon": self.ds.lon},
        )

        weights = get_weights(self.ds, axis=["Y"])
        result = correlation(
            self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["Y"]
        )

        assert_allclose(expected, result)


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

    def test_raises_error_with_invalid_axis_arg(self):
        with pytest.raises(ValueError):
            rmse(self.ds.ts_model, self.ds.ts_obs, "ts", axis=["T"])  # type: ignore

        with pytest.raises(ValueError):
            rmse(self.ds.ts_model, self.ds.ts_obs, "ts", axis="T")  # type: ignore

    def test_returns_weighted_rmse_on_x_y_axes(self):
        expected = xr.DataArray(
            data=np.array([0.13976063, np.nan, 0.07071068], dtype="float64"),
            coords={"time": self.ds.time},
        )

        weights = get_weights(self.ds, axis=["X", "Y"])
        result = rmse(
            self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["X", "Y"]
        )

        assert_allclose(expected, result)

    def test_returns_weighted_rmse_on_x_axis(self):
        expected = xr.DataArray(
            data=np.array(
                [
                    [0.1767767, 0.08838835],
                    [np.nan, 0.07071068],
                    [0.07071068, 0.07071068],
                ],
                dtype="float64",
            ),
            coords={"time": self.ds.time, "lat": self.ds.lat},
        )

        weights = get_weights(self.ds, axis=["X"])
        result = rmse(self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["X"])

        assert_allclose(expected, result)

    def test_returns_weighted_rmse_on_y_axis(self):
        expected = xr.DataArray(
            data=np.array(
                [
                    [0.05302897, 0.19040483],
                    [np.nan, 0.14143213],
                    [0.07070529, 0.07071606],
                ],
                dtype="float64",
            ),
            coords={"time": self.ds.time, "lon": self.ds.lon},
        )

        weights = get_weights(self.ds, axis=["Y"])
        result = rmse(self.ds.ts_model, self.ds.ts_obs, weights=weights, axis=["Y"])

        assert_allclose(expected, result)
