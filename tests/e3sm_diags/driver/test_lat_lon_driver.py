import cftime
import numpy as np
import pytest
import xarray as xr

from e3sm_diags.driver.lat_lon_driver import _get_ref_climo_dataset
from e3sm_diags.driver.utils.dataset_xr import Dataset
from tests.e3sm_diags.driver.utils.test_dataset_xr import (
    _create_parameter_object,
    spatial_bounds,
    spatial_coords,
)


class TestGetReferenceClimoDataset:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        # Create temporary directory to save files.
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        # Set up climatology dataset and save to a temp file.
        self.ds_climo = xr.Dataset(
            coords={
                **spatial_coords,
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 1, 12, 0, 0, 0, has_year_zero=False
                            )
                        ],
                        dtype="object",
                    ),
                    attrs={
                        "axis": "T",
                        "long_name": "time",
                        "standard_name": "time",
                        "bounds": "time_bnds",
                    },
                ),
            },
            data_vars={
                **spatial_bounds,
                "ts": xr.DataArray(
                    name="ts",
                    data=np.array(
                        [
                            [[1.0, 1.0], [1.0, 1.0]],
                        ]
                    ),
                    dims=["time", "lat", "lon"],
                ),
            },
        )
        self.ds_climo.time.encoding = {"units": "days since 2000-01-01"}

        # Set up time series dataset and save to a temp file.
        self.ds_ts = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 1, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 2, 1, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 3, 1, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2001, 1, 1, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype="object",
                    ),
                    attrs={
                        "axis": "T",
                        "long_name": "time",
                        "standard_name": "time",
                        "bounds": "time_bnds",
                    },
                ),
            },
            data_vars={
                "time_bnds": xr.DataArray(
                    name="time_bnds",
                    data=np.array(
                        [
                            [
                                cftime.DatetimeGregorian(
                                    2000, 1, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                                cftime.DatetimeGregorian(
                                    2000, 2, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                            ],
                            [
                                cftime.DatetimeGregorian(
                                    2000, 2, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                                cftime.DatetimeGregorian(
                                    2000, 3, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                            ],
                            [
                                cftime.DatetimeGregorian(
                                    2000, 3, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                                cftime.DatetimeGregorian(
                                    2000, 4, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                            ],
                            [
                                cftime.DatetimeGregorian(
                                    2001, 1, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                                cftime.DatetimeGregorian(
                                    2001, 2, 1, 0, 0, 0, 0, has_year_zero=False
                                ),
                            ],
                        ],
                        dtype=object,
                    ),
                    dims=["time", "bnds"],
                ),
                "ts": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                                [[1.0, 1.0], [1.0, 1.0]],
                                [[1.0, 1.0], [1.0, 1.0]],
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                    )
                ),
            },
        )
        self.ds_ts.time.encoding = {"units": "days since 2000-01-01"}

    def test_raises_error_if_dataset_data_type_is_not_ref(self):
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "test.nc"
        ds = Dataset(parameter, data_type="test")

        with pytest.raises(RuntimeError):
            _get_ref_climo_dataset(ds, "ts", "ANN")

    def test_returns_reference_climo_dataset_from_file(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "ref_file.nc"

        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")
        result = _get_ref_climo_dataset(ds, "ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_None_if_climo_dataset_not_found(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "ref_file.nc"
        ds = Dataset(parameter, data_type="ref")

        result = _get_ref_climo_dataset(ds, "ts", "ANN")

        assert result is None
