from collections import OrderedDict
from typing import Literal

import cftime
import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.acme_new import DERIVED_VARIABLES
from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter


def _create_parameter_object(
    dataset_type: Literal["ref", "test"],
    data_type: Literal["climo", "time_series"],
    data_path: str,
    start_yr: str,
    end_yr: str,
):
    parameter = CoreParameter()

    if dataset_type == "ref":
        if data_type == "time_series":
            parameter.ref_timeseries_input = True
        else:
            parameter.ref_timeseries_input = False

        parameter.reference_data_path = data_path
        parameter.ref_start_yr = start_yr  # type: ignore
        parameter.ref_end_yr = end_yr  # type: ignore
    elif dataset_type == "test":
        if data_type == "time_series":
            parameter.test_timeseries_input = True
        else:
            parameter.test_timeseries_input = False

        parameter.test_data_path = data_path
        parameter.test_start_yr = start_yr  # type: ignore
        parameter.test_end_yr = end_yr  # type: ignore

    return parameter


class TestInit:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

    def test_sets_attrs_if_type_attr_is_ref(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        assert ds.root_path == parameter.reference_data_path
        assert ds.start_yr == parameter.ref_start_yr
        assert ds.end_yr == parameter.ref_end_yr

    def test_sets_attrs_if_type_attr_is_test(self):
        parameter = _create_parameter_object(
            "test", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="test")

        assert ds.root_path == parameter.test_data_path
        assert ds.start_yr == parameter.test_start_yr
        assert ds.end_yr == parameter.test_end_yr

    def test_raises_error_if_type_attr_is_invalid(self):
        parameter = CoreParameter()

        with pytest.raises(ValueError):
            Dataset(parameter, type="invalid")  # type: ignore

    def test_sets_start_yr_and_end_yr_for_area_mean_time_series_set(self):
        parameter = AreaMeanTimeSeriesParameter()
        parameter.sets[0] = "area_mean_time_series"
        parameter.start_yr = "2000"
        parameter.end_yr = "2001"

        ds = Dataset(parameter, type="ref")

        assert ds.start_yr == parameter.start_yr
        assert ds.end_yr == parameter.end_yr

    def test_sets_sub_monthly_if_diurnal_cycle_or_arms_diags_set(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.sets[0] = "diurnal_cycle"

        ds = Dataset(parameter, type="ref")

        assert ds.is_sub_monthly

        parameter.sets[0] = "arm_diags"
        ds2 = Dataset(parameter, type="ref")

        assert ds2.is_sub_monthly

    def test_sets_derived_vars_map(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        assert ds.derived_vars_map == DERIVED_VARIABLES

    def test_sets_drived_vars_map_with_existing_entry(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.derived_variables = {
            "PRECT": OrderedDict([(("some_var",), lambda some_var: some_var)])
        }

        ds = Dataset(parameter, type="ref")

        # The expected `derived_vars_map` result.
        expected = DERIVED_VARIABLES.copy()
        expected["PRECT"] = OrderedDict(
            **parameter.derived_variables["PRECT"], **expected["PRECT"]
        )

        assert ds.derived_vars_map == expected

    def test_sets_drived_vars_map_with_new_entry(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.derived_variables = {
            "NEW_DERIVED_VAR": OrderedDict([(("some_var",), lambda some_var: some_var)])
        }

        ds = Dataset(parameter, type="ref")

        # The expected `derived_vars_map` result.
        expected = DERIVED_VARIABLES.copy()
        expected["NEW_DERIVED_VAR"] = parameter.derived_variables["NEW_DERIVED_VAR"]

        assert ds.derived_vars_map == expected


class TestDataSetProperties:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

    def test_property_is_timeseries_returns_true_and_is_climo_returns_false_for_ref(
        self,
    ):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_true_and_is_climo_returns_false_for_test(
        self,
    ):
        parameter = _create_parameter_object(
            "test", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="test")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_test(
        self,
    ):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        assert not ds.is_time_series
        assert ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_ref(
        self,
    ):
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        ds = Dataset(parameter, type="test")

        assert not ds.is_time_series
        assert ds.is_climo


class TestGetClimoDataset:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        # Create temporary directory to save files.
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        # Set up climatology dataset and save to a temp file.
        # TODO: Update this to an actual climatology dataset structure
        self.ds_climo = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
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
                "ts": xr.DataArray(
                    name="ts",
                    data=np.array(
                        [
                            [[1.0, 1.0], [1.0, 1.0]],
                        ]
                    ),
                    dims=["time", "lat", "lon"],
                )
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

    def test_raises_error_if_var_arg_is_not_valid(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var=1, season="ANN")  # type: ignore

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var="", season="ANN")

    def test_raises_error_if_season_arg_is_not_valid(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var="PRECT", season="invalid_season")  # type: ignore

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var="PRECT", season=1)  # type: ignore

    def test_returns_climo_dataset_using_ref_file_variable(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "ref_file.nc"

        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_test_file_variable(self):
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_file = "test_file.nc"

        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.test_file}")

        ds = Dataset(parameter, type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_ref_file_variable_test_name_and_season(self):
        # Example: {test_data_path}/{test_name}_{season}.nc
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_name = "historical_H1"
        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_name}_ANN.nc")

        ds = Dataset(parameter, type="ref")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_test_file_variable_test_name_and_season(self):
        # Example: {test_data_path}/{test_name}_{season}.nc
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_name = "historical_H1"
        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.test_name}_ANN.nc")

        ds = Dataset(parameter, type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_test_file_variable_ref_name_and_season_nested_pattern_1(
        self,
    ):
        # Example: {test_data_path}/{test_name}/{test_name}_{season}.nc
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_name = "historical_H1"

        nested_root_path = self.data_path / parameter.test_name
        nested_root_path.mkdir()

        self.ds_climo.to_netcdf(f"{nested_root_path}/{parameter.test_name}_ANN.nc")

        ds = Dataset(parameter, type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_test_file_variable_ref_name_and_season_nested_pattern_2(
        self,
    ):
        # Example: {test_data_path}/{test_name}/{test_name}_<N CHARACTERS>_{season}.nc
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_name = "historical_H1"

        nested_root_path = self.data_path / parameter.test_name
        nested_root_path.mkdir()

        self.ds_climo.to_netcdf(
            f"{nested_root_path}/{parameter.test_name}_some_other_info_ANN.nc"
        )

        ds = Dataset(parameter, type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_with_derived_variable(self):
        # We will derive the "PRECT" variable using the "pr" variable.
        ds_pr = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "pr": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                        attrs={"units": "mm/s"},
                    )
                ),
            },
        )

        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "pr_200001_200112.nc"
        ds_pr.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        result = ds.get_climo_dataset("PRECT", season="ANN")
        expected = ds_pr.copy()
        expected = expected.squeeze(dim="time")
        expected["PRECT"] = expected["pr"] * 3600 * 24
        expected["PRECT"].attrs["units"] = "mm/day"

        assert result.identical(expected)

    def test_returns_climo_dataset_using_derived_var_directly_from_dataset(self):
        ds_precst = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "PRECST": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                        attrs={"units": "mm/s"},
                    )
                ),
            },
        )

        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "pr_200001_200112.nc"
        ds_precst.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        result = ds.get_climo_dataset("PRECST", season="ANN")
        expected = ds_precst.squeeze(dim="time")

        assert result.identical(expected)

    def test_returns_climo_dataset_using_source_variable_with_wildcard(self):
        ds_precst = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "bc_a?DDF": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                    )
                ),
                "bc_c?DDF": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                    )
                ),
            },
        )

        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "var_200001_200112.nc"
        ds_precst.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        result = ds.get_climo_dataset("bc_DDF", season="ANN")
        expected = ds_precst.squeeze(dim="time")
        expected["bc_DDF"] = expected["bc_a?DDF"] + expected["bc_c?DDF"]

        assert result.identical(expected)

    def test_returns_climo_dataset_using_climo_of_time_series_files(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.ref_timeseries_input = True
        parameter.ref_file = "ts_200001_201112.nc"

        self.ds_ts.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        result = ds.get_climo_dataset("ts", "ANN")
        # Since the data is not sub-monthly, the first time coord (2001-01-01)
        # is dropped when subsetting with the middle of the month (2000-01-15).
        expected = self.ds_ts.isel(time=slice(1, 4))
        expected["ts"] = xr.DataArray(
            name="ts", data=np.array([[1.0, 1.0], [1.0, 1.0]]), dims=["lat", "lon"]
        )

        assert result.identical(expected)

    def test_raises_error_if_no_filepath_found_for_variable(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        parameter.ref_timeseries_input = False

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("some_var", "ANN")

    def test_raises_error_if_var_not_in_dataset_or_derived_var_map(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        parameter.ref_timeseries_input = False
        parameter.ref_file = "ts_200001_201112.nc"

        self.ds_ts.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("some_var", "ANN")

    def test_raises_error_if_dataset_has_no_matching_source_variables_to_derive_variable(
        self,
    ):
        # In this test, we don't create a dataset and write it out to `.nc`.
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "pr_200001_200112.nc"

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("PRECT", season="ANN")

    def test_raises_error_if_no_datasets_found_to_derive_variable(self):
        ds_precst = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "invalid": xr.DataArray(
                    xr.DataArray(
                        data=np.array(
                            [
                                [[1.0, 1.0], [1.0, 1.0]],
                            ]
                        ),
                        dims=["time", "lat", "lon"],
                        attrs={"units": "mm/s"},
                    )
                ),
            },
        )

        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "pr_200001_200112.nc"
        ds_precst.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("PRECST", season="ANN")


class TestGetTimeSeriesDataset:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        # Set up time series dataset and save to a temp file.
        self.ts_path = f"{self.data_path}/ts_200001_200112.nc"
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
                        dtype=object,
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

    def test_raises_error_if_data_is_not_time_series(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.ref_timeseries_input = False

        ds = Dataset(parameter, type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset(var="ts")

    def test_raises_error_if_var_arg_is_not_valid(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        # Not a string
        with pytest.raises(ValueError):
            ds.get_time_series_dataset(var=1)  # type: ignore

        # An empty string
        with pytest.raises(ValueError):
            ds.get_time_series_dataset(var="")

    def test_returns_time_series_dataset_using_file(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        result = ds.get_time_series_dataset("ts")

        # Since the data is not sub-monthly, the first time coord (2001-01-01)
        # is dropped when subsetting with the middle of the month (2000-01-15).
        expected = self.ds_ts.isel(time=slice(1, 4))

        assert result.identical(expected)

    def test_returns_time_series_dataset_using_sub_monthly_sets(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        self.ds_ts.to_netcdf(f"{self.data_path}/ts_200001_200112.nc")
        # "arm_diags" includes the the regions parameter in the filename
        self.ds_ts.to_netcdf(f"{self.data_path}/ts_global_200001_200112.nc")

        for set in ["diurnal_cycle", "arm_diags"]:
            parameter.sets[0] = set

            ds = Dataset(parameter, type="ref")

            result = ds.get_time_series_dataset("ts")
            expected = self.ds_ts.copy()

            assert result.identical(expected)

    def test_returns_time_series_dataset_using_derived_var(self):
        # We will derive the "PRECT" variable using the "pr" variable.
        ds_pr = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 2, 15, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 3, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2001, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "pr": xr.DataArray(
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
                        attrs={"units": "mm/s"},
                    )
                ),
            },
        )
        ds_pr.to_netcdf(f"{self.data_path}/pr_200001_200112.nc")

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        result = ds.get_time_series_dataset("PRECT")
        expected = ds_pr.copy()
        expected["PRECT"] = expected["pr"] * 3600 * 24
        expected["PRECT"].attrs["units"] = "mm/day"

        assert result.identical(expected)

    def test_returns_time_series_dataset_using_derived_var_directly_from_dataset(self):
        # We will derive the "PRECT" variable using the "pr" variable.
        ds_precst = xr.Dataset(
            coords={
                "lat": [-90, 90],
                "lon": [0, 180],
                "time": xr.DataArray(
                    dims="time",
                    data=np.array(
                        [
                            cftime.DatetimeGregorian(
                                2000, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 2, 15, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2000, 3, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                            cftime.DatetimeGregorian(
                                2001, 1, 16, 12, 0, 0, 0, has_year_zero=False
                            ),
                        ],
                        dtype=object,
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
                "PRECST": xr.DataArray(
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
                        attrs={"units": "mm/s"},
                    )
                ),
            },
        )
        ds_precst.to_netcdf(f"{self.data_path}/PRECST_200001_200112.nc")

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        result = ds.get_time_series_dataset("PRECST")
        expected = ds_precst.copy()

        assert result.identical(expected)

    def test_raises_error_if_no_datasets_found_to_derive_variable(self):
        # In this test, we don't create a dataset and write it out to `.nc`.
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("PRECT")

    def test_returns_time_series_dataset_with_centered_time_if_single_point(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.sets[0] = "diurnal_cycle"

        ds = Dataset(parameter, type="ref")

        result = ds.get_time_series_dataset("ts", single_point=True)
        expected = self.ds_ts.copy()
        expected["time"].data[:] = np.array(
            [
                cftime.DatetimeGregorian(2000, 1, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 2, 15, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 3, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2001, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            ],
            dtype=object,
        )

        assert result.identical(expected)

    def test_returns_time_series_dataset_using_file_with_ref_name_prepended(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.ref_name = "historical_H1"

        ref_data_path = self.data_path / parameter.ref_name
        ref_data_path.mkdir()
        self.ds_ts.to_netcdf(f"{ref_data_path}/ts_200001_200112.nc")

        ds = Dataset(parameter, type="ref")

        result = ds.get_time_series_dataset("ts")
        # Since the data is not sub-monthly, the first time coord (2001-01-01)
        # is dropped when subsetting with the middle of the month (2000-01-15).
        expected = self.ds_ts.isel(time=slice(1, 4))

        assert result.identical(expected)

    def test_raises_error_if_time_series_dataset_could_not_be_found(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("invalid_var")

    def test_raises_error_if_multiple_time_series_datasets_found_for_single_var(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        self.ds_ts.to_netcdf(f"{self.data_path}/ts_199901_200012.nc")
        ds = Dataset(parameter, type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("ts")

    def test_raises_error_when_time_slicing_if_start_year_less_than_var_start_year(
        self,
    ):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "1999", "2001"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset("ts")

    def test_raises_error_when_time_slicing_if_end_year_greater_than_var_end_year(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2002"
        )

        ds = Dataset(parameter, type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset("ts")


class Test_GetLandSeaMask:
    def test_returns_land_sea_mask_if_matching_vars_in_dataset(self):
        assert 0

    def test_returns_default_land_sea_mask_if_one_or_no_matching_vars_in_dataset(self):
        assert 0

    def test_drops_time_dimension_if_time_dimension_in_land_sea_mask_dataset(self):
        assert 0
