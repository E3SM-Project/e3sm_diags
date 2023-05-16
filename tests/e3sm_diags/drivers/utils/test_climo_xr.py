import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver.utils.climo_xr import climo


class TestClimo:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        # Create temporary directory to save files.
        dir = tmp_path / "input_data"
        dir.mkdir()

        ds = xr.Dataset(
            data_vars={
                "ts": xr.DataArray(
                    data=np.array(
                        [[[2.0]], [[1.0]], [[1.0]], [[1.0]], [[2.0]]], dtype="float64"
                    ),
                    dims=["time", "lat", "lon"],
                    attrs={"test_attr": "test"},
                ),
                "time_bnds": xr.DataArray(
                    name="time_bnds",
                    data=np.array(
                        [
                            [
                                "2000-01-01T00:00:00.000000000",
                                "2000-02-01T00:00:00.000000000",
                            ],
                            [
                                "2000-03-01T00:00:00.000000000",
                                "2000-04-01T00:00:00.000000000",
                            ],
                            [
                                "2000-06-01T00:00:00.000000000",
                                "2000-07-01T00:00:00.000000000",
                            ],
                            [
                                "2000-09-01T00:00:00.000000000",
                                "2000-10-01T00:00:00.000000000",
                            ],
                            [
                                "2001-02-01T00:00:00.000000000",
                                "2001-03-01T00:00:00.000000000",
                            ],
                        ],
                        dtype="datetime64[ns]",
                    ),
                    dims=["time", "bnds"],
                    attrs={"xcdat_bounds": "True"},
                ),
            },
            coords={
                "lat": xr.DataArray(
                    data=np.array([-90]),
                    dims=["lat"],
                    attrs={
                        "axis": "Y",
                        "long_name": "latitude",
                        "standard_name": "latitude",
                    },
                ),
                "lon": xr.DataArray(
                    data=np.array([0]),
                    dims=["lon"],
                    attrs={
                        "axis": "X",
                        "long_name": "longitude",
                        "standard_name": "longitude",
                    },
                ),
                "time": xr.DataArray(
                    data=np.array(
                        [
                            "2000-01-16T12:00:00.000000000",
                            "2000-03-16T12:00:00.000000000",
                            "2000-06-16T00:00:00.000000000",
                            "2000-09-16T00:00:00.000000000",
                            "2001-02-15T12:00:00.000000000",
                        ],
                        dtype="datetime64[ns]",
                    ),
                    dims=["time"],
                    attrs={
                        "axis": "T",
                        "long_name": "time",
                        "standard_name": "time",
                        "bounds": "time_bnds",
                    },
                ),
            },
        )
        ds.time.encoding = {"units": "days since 2000-01-01"}

        # Write the dataset to an `.nc` file and set the DataArray encoding
        # attribute to mimic a real-world dataset.
        filepath = f"{dir}/file.nc"
        ds.to_netcdf(filepath)
        ds.ts.encoding["source"] = filepath

        self.ts = ds.ts.copy()

    def test_returns_annual_cycle_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "ANN")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.42309607]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_DJF_season_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "DJF")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[2.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_MAM_season_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "MAM")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_JJA_season_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "JJA")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_SON_season_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "SON")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_jan_climatology(self):
        ts = self.ts.copy()

        result = climo(ts, "01")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[2.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs

    def test_returns_climatology_for_derived_variable(self):
        ts = self.ts.copy()

        # Delete the source of this variable to mimic a "derived" variable,
        # which is a variable created using other variables in the dataset.
        del ts.encoding["source"]

        result = climo(ts, "01")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[2.0]]),
            coords={
                "lat": ts.lat,
                "lon": ts.lon,
            },
            dims=["lat", "lon"],
            attrs={"test_attr": "test"},
        )

        # Check DataArray values and attributes align
        xr.testing.assert_allclose(result, expected)
        assert result.attrs == expected.attrs

        for coord in result.coords:
            assert result[coord].attrs == expected[coord].attrs
