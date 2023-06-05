from pathlib import Path

import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver.utils.climo_xcdat import climo


class TestClimo:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path: Path):
        # Create temporary directory to save files.
        dir = tmp_path / "input_data"
        dir.mkdir()

        self.ds = xr.Dataset(
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
        self.ds.time.encoding = {
            "units": "years since 2000-01-01",
            "calendar": "standard",
        }

    def test_returns_annual_cycle_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "ANN")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.4]]),
            dims=["lat", "lon"],
            coords={"lat": ds.lat, "lon": ds.lon},
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "average",
                "freq": "month",
                "weighted": "True",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)

    def test_returns_DJF_season_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "DJF")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[2.0]]),
            coords={
                "lat": ds.lat,
                "lon": ds.lon,
            },
            dims=["lat", "lon"],
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "climatology",
                "freq": "season",
                "weighted": "True",
                "dec_mode": "DJF",
                "drop_incomplete_djf": "False",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)

    def test_returns_MAM_season_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "MAM")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ds.lat,
                "lon": ds.lon,
            },
            dims=["lat", "lon"],
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "climatology",
                "freq": "season",
                "dec_mode": "DJF",
                "drop_incomplete_djf": "False",
                "weighted": "True",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)

    def test_returns_JJA_season_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "JJA")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ds.lat,
                "lon": ds.lon,
            },
            dims=["lat", "lon"],
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "climatology",
                "freq": "season",
                "dec_mode": "DJF",
                "drop_incomplete_djf": "False",
                "weighted": "True",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)

    def test_returns_SON_season_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "SON")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[1.0]]),
            coords={
                "lat": ds.lat,
                "lon": ds.lon,
            },
            dims=["lat", "lon"],
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "climatology",
                "freq": "season",
                "dec_mode": "DJF",
                "drop_incomplete_djf": "False",
                "weighted": "True",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)

    def test_returns_jan_climatology(self):
        ds = self.ds.copy()

        result = climo(ds, "ts", "01")
        expected = xr.DataArray(
            name="ts",
            data=np.array([[2.0]]),
            coords={
                "lat": ds.lat,
                "lon": ds.lon,
            },
            dims=["lat", "lon"],
            attrs={
                "test_attr": "test",
                "operation": "temporal_avg",
                "mode": "climatology",
                "freq": "month",
                "weighted": "True",
            },
        )

        # Check DataArray values and attributes align
        assert result.identical(expected)
