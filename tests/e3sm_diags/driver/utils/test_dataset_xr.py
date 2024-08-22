import copy
import logging
from collections import OrderedDict
from typing import Literal

import cftime
import numpy as np
import pytest
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
from e3sm_diags.driver import LAND_OCEAN_MASK_PATH
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter

# Reusable spatial coords dictionary for composing an xr.Dataest.
spatial_coords = {
    "lat": xr.DataArray(
        dims="lat",
        data=np.array([-90.0, 90]),
        attrs={
            "axis": "Y",
            "long_name": "latitude",
            "standard_name": "latitude",
            "bounds": "lat_bnds",
        },
    ),
    "lon": xr.DataArray(
        dims="lon",
        data=np.array([0.0, 180]),
        attrs={
            "axis": "X",
            "long_name": "longitude",
            "standard_name": "longitude",
            "bounds": "lon_bnds",
        },
    ),
}

# Reusable spatial bounds dictionary for composing an xr.Dataest.
spatial_bounds = {
    "lat_bnds": xr.DataArray(
        name="lat_bnds",
        data=[[-90.0, 0.0], [0.0, 90.0]],
        dims=["lat", "bnds"],
    ),
    "lon_bnds": xr.DataArray(
        name="lat_bnds",
        data=[[-90.0, 90.0], [90, 270]],
        dims=["lon", "bnds"],
    ),
}


def _create_parameter_object(
    dataset_type: Literal["ref", "test"],
    data_type: Literal["climo", "time_series"],
    data_path: str,
    start_yr: str,
    end_yr: str,
):
    # NOTE: Make sure to create deep copies to avoid references in memory to
    # the same object.
    parameter = copy.deepcopy(CoreParameter())

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

        ds = Dataset(parameter, data_type="ref")

        assert ds.root_path == parameter.reference_data_path
        assert ds.start_yr == parameter.ref_start_yr
        assert ds.end_yr == parameter.ref_end_yr

    def test_sets_attrs_if_type_attr_is_test(self):
        parameter = _create_parameter_object(
            "test", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="test")

        assert ds.root_path == parameter.test_data_path
        assert ds.start_yr == parameter.test_start_yr
        assert ds.end_yr == parameter.test_end_yr

    def test_raises_error_if_type_attr_is_invalid(self):
        parameter = CoreParameter()

        with pytest.raises(ValueError):
            Dataset(parameter, data_type="invalid")  # type: ignore

    def test_sets_start_yr_and_end_yr_for_area_mean_time_series_set(self):
        parameter = AreaMeanTimeSeriesParameter()
        parameter.sets = ["area_mean_time_series"]
        parameter.start_yr = "2000"
        parameter.end_yr = "2001"

        ds = Dataset(parameter, data_type="ref")

        assert ds.start_yr == parameter.start_yr
        assert ds.end_yr == parameter.end_yr

    def test_sets_sub_monthly_if_diurnal_cycle_or_arms_diags_set(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.sets = ["diurnal_cycle"]

        ds = Dataset(parameter, data_type="ref")

        assert ds.is_sub_monthly

        parameter.sets[0] = "arm_diags"
        ds2 = Dataset(parameter, data_type="ref")

        assert ds2.is_sub_monthly

    def test_sets_derived_vars_map(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        assert ds.derived_vars_map == DERIVED_VARIABLES

    def test_sets_drived_vars_map_with_existing_entry(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.derived_variables = {
            "PRECT": OrderedDict([(("some_var",), lambda some_var: some_var)])
        }

        ds = Dataset(parameter, data_type="ref")

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

        ds = Dataset(parameter, data_type="ref")

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

        ds = Dataset(parameter, data_type="ref")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_true_and_is_climo_returns_false_for_test(
        self,
    ):
        parameter = _create_parameter_object(
            "test", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="test")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_test(
        self,
    ):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        assert not ds.is_time_series
        assert ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_ref(
        self,
    ):
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        ds = Dataset(parameter, data_type="test")

        assert not ds.is_time_series
        assert ds.is_climo


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
            ds.get_ref_climo_dataset("ts", "ANN", self.ds_climo.copy())

    def test_returns_reference_climo_dataset_from_file(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "ref_file.nc"

        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")
        result = ds.get_ref_climo_dataset("ts", "ANN", self.ds_climo.copy())
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)
        assert not ds.model_only

    def test_returns_test_dataset_as_default_value_if_climo_dataset_not_found(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_file = "ref_file.nc"
        ds = Dataset(parameter, data_type="ref")

        ds_test = self.ds_climo.copy()
        result = ds.get_ref_climo_dataset("ts", "ANN", ds_test)

        assert result.identical(ds_test)
        assert ds.model_only


class TestGetClimoDataset:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        # Create temporary directory to save files.
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        self.spatial_coords = {
            "lat": xr.DataArray(
                dims="lat",
                data=np.array([-90.0, 90]),
                attrs={
                    "axis": "Y",
                    "long_name": "latitude",
                    "standard_name": "latitude",
                    "bounds": "lat_bnds",
                },
            ),
            "lon": xr.DataArray(
                dims="lon",
                data=np.array([0.0, 180]),
                attrs={
                    "axis": "X",
                    "long_name": "longitude",
                    "standard_name": "longitude",
                    "bounds": "lon_bnds",
                },
            ),
        }
        self.spatial_bounds = {
            "lat_bnds": xr.DataArray(
                name="lat_bnds",
                data=[[-90.0, 0.0], [0.0, 90.0]],
                dims=["lat", "bnds"],
            ),
            "lon_bnds": xr.DataArray(
                name="lat_bnds",
                data=[[-90.0, 90.0], [90, 270]],
                dims=["lon", "bnds"],
            ),
        }

        # Set up climatology dataset and save to a temp file.
        # TODO: Update this to an actual climatology dataset structure
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
                **spatial_coords,
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

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var=1, season="ANN")  # type: ignore

        with pytest.raises(ValueError):
            ds.get_climo_dataset(var="", season="ANN")

    def test_raises_error_if_season_arg_is_not_valid(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

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

        ds = Dataset(parameter, data_type="ref")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_test_file_variable(self):
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_file = "test_file.nc"

        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.test_file}")

        ds = Dataset(parameter, data_type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_ref_file_variable_test_name_and_season(self):
        # Example: {test_data_path}/{test_name}_{season}.nc
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )
        parameter.ref_name = "historical_H1"
        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_name}_ANN.nc")

        ds = Dataset(parameter, data_type="ref")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_test_file_variable_test_name_and_season(self):
        # Example: {test_data_path}/{test_name}_{season}.nc
        parameter = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2001"
        )
        parameter.test_name = "historical_H1"
        self.ds_climo.to_netcdf(f"{self.data_path}/{parameter.test_name}_ANN.nc")

        ds = Dataset(parameter, data_type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

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

        ds = Dataset(parameter, data_type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

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

        ds = Dataset(parameter, data_type="test")
        result = ds.get_climo_dataset("ts", "ANN")
        expected = self.ds_climo.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_with_derived_variable(self):
        # We will derive the "PRECT" variable using the "pr" variable.
        ds_pr = xr.Dataset(
            coords={
                **spatial_coords,
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
                **spatial_bounds,
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

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_climo_dataset("PRECT", season="ANN")
        expected = ds_pr.copy()
        expected = expected.squeeze(dim="time").drop_vars("time")
        expected["PRECT"] = expected["pr"] * 3600 * 24
        expected["PRECT"].attrs["units"] = "mm/day"
        expected = expected.drop_vars("pr")

        xr.testing.assert_identical(result, expected)

    @pytest.mark.xfail
    def test_returns_climo_dataset_using_derived_var_directly_from_dataset_and_replaces_scalar_time_var(
        self,
    ):
        # FIXME: This test needs to cover `except` block in `_open_dataset()`.
        # The issue is that we can't create a dummy dataset with an incorrect
        # time scalar variable using Xarray because it just throws the error
        # below. We might need to use another library like netCDF4 to create
        # a dummy dataset.
        ds_src = xr.Dataset(
            coords={
                **spatial_coords,
            },
            data_vars={
                **spatial_bounds,
                "time": xr.DataArray(
                    dims="time",
                    data=0,
                ),
                "SOURCE_VAR": xr.DataArray(
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
        ds_src.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_climo_dataset("SO", season="ANN")
        expected = ds_src.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_derived_var_directly_from_dataset(self):
        ds_src = xr.Dataset(
            coords={
                **spatial_coords,
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
                **spatial_bounds,
                "SOURCE_VAR": xr.DataArray(
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
        ds_src.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_climo_dataset("SOURCE_VAR", season="ANN")
        expected = ds_src.squeeze(dim="time").drop_vars("time")

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_source_variable_with_wildcard(self):
        ds_precst = xr.Dataset(
            coords={
                **spatial_coords,
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
                **spatial_bounds,
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

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_climo_dataset("bc_DDF", season="ANN")
        expected = ds_precst.squeeze(dim="time").drop_vars("time")
        expected["bc_DDF"] = expected["bc_a?DDF"] + expected["bc_c?DDF"]
        expected = expected.drop_vars(["bc_a?DDF", "bc_c?DDF"])

        xr.testing.assert_identical(result, expected)

    def test_returns_climo_dataset_using_climo_of_time_series_files(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.ref_timeseries_input = True
        parameter.ref_file = "ts_200001_201112.nc"

        self.ds_ts.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_climo_dataset("ts", "ANN")

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
        expected["ts"] = xr.DataArray(
            name="ts", data=np.array([[1.0, 1.0], [1.0, 1.0]]), dims=["lat", "lon"]
        )
        # Set all of the correct attributes.
        expected = expected.assign(**spatial_coords, **spatial_bounds)  # type: ignore
        expected["lat"].attrs["units"] = "degrees_north"
        expected["lat_bnds"].attrs["xcdat_bounds"] = "True"
        expected["lon_bnds"].attrs["xcdat_bounds"] = "True"

        xr.testing.assert_identical(result, expected)

    def test_raises_error_if_no_filepath_found_for_variable(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        parameter.ref_timeseries_input = False

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("some_var", "ANN")

    def test_raises_error_if_var_not_in_dataset_or_derived_var_map(self):
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2001"
        )

        parameter.ref_timeseries_input = False
        parameter.ref_file = "ts_200001_201112.nc"

        self.ds_ts.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")

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

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(IOError):
            ds.get_climo_dataset("PRECT", season="ANN")

    def test_raises_error_if_no_datasets_found_to_derive_variable(self):
        ds_precst = xr.Dataset(
            coords={
                **spatial_coords,
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
                **spatial_bounds,
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

        ds = Dataset(parameter, data_type="ref")

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
                **spatial_coords,
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
                **spatial_bounds,
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

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset(var="ts")

    def test_raises_error_if_var_arg_is_not_valid(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

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

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_time_series_dataset("ts")

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

        xr.testing.assert_identical(result, expected)

    def test_returns_time_series_dataset_using_sub_monthly_sets(self):
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        self.ds_ts.to_netcdf(f"{self.data_path}/ts_200001_200112.nc")
        # "arm_diags" includes the the regions parameter in the filename
        self.ds_ts.to_netcdf(f"{self.data_path}/ts_global_200001_200112.nc")

        for set in ["diurnal_cycle", "arm_diags"]:
            parameter.sets[0] = set

            ds = Dataset(parameter, data_type="ref")

            result = ds.get_time_series_dataset("ts")
            expected = self.ds_ts.copy()

            xr.testing.assert_identical(result, expected)

    def test_returns_time_series_dataset_using_derived_var(self):
        # We will derive the "PRECT" variable using the "pr" variable.
        ds_pr = self.ds_ts.copy()
        ds_pr["pr"] = xr.DataArray(
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
        ds_pr = ds_pr.drop_vars("ts")
        ds_pr.to_netcdf(f"{self.data_path}/pr_200001_200112.nc")

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_time_series_dataset("PRECT")
        expected = ds_pr.copy()
        expected["time"].data[:] = np.array(
            [
                cftime.DatetimeGregorian(2000, 1, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 2, 15, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 3, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2001, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            ],
            dtype=object,
        )
        expected["PRECT"] = expected["pr"] * 3600 * 24
        expected["PRECT"].attrs["units"] = "mm/day"

        xr.testing.assert_identical(result, expected)

    def test_returns_time_series_dataset_using_derived_var_directly_from_dataset(self):
        ds_precst = self.ds_ts.copy()
        ds_precst["PRECST"] = xr.DataArray(
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
        ds_precst = ds_precst.drop_vars("ts")
        ds_precst.to_netcdf(f"{self.data_path}/PRECST_200001_200112.nc")

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_time_series_dataset("PRECST")
        expected = ds_precst.copy()
        expected = ds_precst.copy()
        expected["PRECST"].attrs["units"] = "mm/s"
        expected["time"].data[:] = np.array(
            [
                cftime.DatetimeGregorian(2000, 1, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 2, 15, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 3, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2001, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            ],
            dtype=object,
        )

        xr.testing.assert_identical(result, expected)

    def test_raises_error_if_no_datasets_found_to_derive_variable(self):
        # In this test, we don't create a dataset and write it out to `.nc`.
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("PRECT")

    def test_returns_time_series_dataset_without_centered_time_if_single_point_data(
        self,
    ):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.sets[0] = "diurnal_cycle"

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_time_series_dataset("ts", single_point=True)
        expected = self.ds_ts.copy()

        xr.testing.assert_identical(result, expected)

    def test_returns_time_series_dataset_with_centered_time_if_non_sub_monthly_data(
        self,
    ):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")
        ds.is_sub_monthly = False

        result = ds.get_time_series_dataset("ts")
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

        xr.testing.assert_identical(result, expected)

    def test_returns_time_series_dataset_using_file_with_ref_name_prepended(self):
        ds_ts = self.ds_ts.copy()
        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        parameter.ref_name = "historical_H1"

        ref_data_path = self.data_path / parameter.ref_name
        ref_data_path.mkdir()
        ds_ts.to_netcdf(f"{ref_data_path}/ts_200001_200112.nc")

        ds = Dataset(parameter, data_type="ref")

        result = ds.get_time_series_dataset("ts")

        expected = ds_ts.copy()
        expected["time"].data[:] = np.array(
            [
                cftime.DatetimeGregorian(2000, 1, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 2, 15, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2000, 3, 16, 12, 0, 0, 0, has_year_zero=False),
                cftime.DatetimeGregorian(2001, 1, 16, 12, 0, 0, 0, has_year_zero=False),
            ],
            dtype=object,
        )

        xr.testing.assert_identical(result, expected)

    def test_raises_error_if_time_series_dataset_could_not_be_found(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("invalid_var")

    def test_raises_error_if_multiple_time_series_datasets_found_for_single_var(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2001"
        )
        self.ds_ts.to_netcdf(f"{self.data_path}/ts_199901_200012.nc")
        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(IOError):
            ds.get_time_series_dataset("ts")

    def test_raises_error_when_time_slicing_if_start_year_less_than_var_start_year(
        self,
    ):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "1999", "2001"
        )

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset("ts")

    def test_raises_error_when_time_slicing_if_end_year_greater_than_var_end_year(self):
        self.ds_ts.to_netcdf(self.ts_path)

        parameter = _create_parameter_object(
            "ref", "time_series", self.data_path, "2000", "2002"
        )

        ds = Dataset(parameter, data_type="ref")

        with pytest.raises(ValueError):
            ds.get_time_series_dataset("ts")


class Test_GetLandSeaMask:
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

    def test_returns_land_sea_mask_if_matching_vars_in_dataset(self):
        ds_climo: xr.Dataset = self.ds_climo.copy()
        ds_climo["LANDFRAC"] = xr.DataArray(
            name="LANDFRAC",
            data=[
                [[1.0, 1.0], [1.0, 1.0]],
            ],
            dims=["time", "lat", "lon"],
        )
        ds_climo["OCNFRAC"] = xr.DataArray(
            name="OCNFRAC",
            data=[
                [[1.0, 1.0], [1.0, 1.0]],
            ],
            dims=["time", "lat", "lon"],
        )

        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2002"
        )
        parameter.ref_file = "ref_file.nc"

        ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")
        result = ds._get_land_sea_mask("ANN")
        expected = ds_climo.copy()
        expected = expected.squeeze(dim="time").drop_vars("time")
        expected = expected.drop_vars("ts")

        xr.testing.assert_identical(result, expected)

    def test_returns_default_land_sea_mask_if_one_or_no_matching_vars_in_dataset(
        self, caplog
    ):
        # Silence logger warning to not pollute test suite.
        caplog.set_level(logging.CRITICAL)

        ds_climo: xr.Dataset = self.ds_climo.copy()
        parameter = _create_parameter_object(
            "ref", "climo", self.data_path, "2000", "2002"
        )
        parameter.ref_file = "ref_file.nc"

        ds_climo.to_netcdf(f"{self.data_path}/{parameter.ref_file}")

        ds = Dataset(parameter, data_type="ref")
        result = ds._get_land_sea_mask("ANN")

        expected = xr.open_dataset(LAND_OCEAN_MASK_PATH)
        expected = expected.squeeze(dim="time").drop_vars(["time", "time_bnds"])

        xr.testing.assert_identical(result, expected)


class TestGetNameAndYearsAttr:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        self.ts_path = f"{self.data_path}/ts_200001_200112.nc"

        # Used for getting climo dataset via `parameter.ref_file`.
        self.ref_file = " ref_file.nc"
        self.test_file = "test_file.nc"

        self.ds_climo = xr.Dataset(attrs={"yrs_averaged": "2000-2002"})

        self.ds_ts = xr.Dataset()
        self.ds_ts.to_netcdf(self.ts_path)

    def test_raises_error_if_season_arg_is_not_passed_for_climo_dataset(self):
        param1 = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2002"
        )
        param1.short_test_name = "short_test_name"

        ds1 = Dataset(param1, data_type="test")

        with pytest.raises(ValueError):
            ds1.get_name_yrs_attr()

    def test_returns_test_name_and_yrs_averaged_attr_with_climo_dataset_using_short_test_name(
        self,
    ):
        # Case 1: name is taken from `parameter.short_test_name`
        param = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2002"
        )
        param.short_test_name = "short_test_name"
        param.test_file = self.test_file

        # Write the climatology dataset out before function call.
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.test_file}")

        ds1 = Dataset(param, data_type="test")
        result = ds1.get_name_yrs_attr("ANN")
        expected = "short_test_name (2000-2002)"

        assert result == expected

    def test_returns_test_name_and_yrs_averaged_attr_with_climo_dataset_using_test_name(
        self,
    ):
        # Case 2: name is taken from `parameter.test_name`
        param = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2002"
        )
        param.test_name = "test_name"

        # Write the climatology dataset out before function call.
        param.test_file = self.test_file
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.test_file}")

        ds2 = Dataset(param, data_type="test")
        result = ds2.get_name_yrs_attr("ANN")
        expected = "test_name (2000-2002)"

        assert result == expected

    def test_returns_only_test_name_attr_if_yrs_averaged_attr_not_found_with_climo_dataset(
        self,
    ):
        param1 = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2002"
        )
        param1.short_test_name = "short_test_name"
        param1.test_file = self.test_file

        # Write the climatology dataset out before function call.
        ds_climo = self.ds_climo.copy()
        del ds_climo.attrs["yrs_averaged"]
        ds_climo.to_netcdf(f"{self.data_path}/{param1.test_file}")

        ds1 = Dataset(param1, data_type="test")
        result = ds1.get_name_yrs_attr("ANN")
        expected = "short_test_name"

        assert result == expected

    def test_returns_only_yrs_averaged_attr_if_test_name_not_set_with_climo_dataset(
        self,
    ):
        param1 = _create_parameter_object(
            "test", "climo", self.data_path, "2000", "2002"
        )
        param1.test_name = ""
        param1.test_file = self.test_file

        # Write the climatology dataset out before function call.
        ds_climo = self.ds_climo.copy()
        ds_climo.to_netcdf(f"{self.data_path}/{param1.test_file}")

        ds1 = Dataset(param1, data_type="test")
        result = ds1.get_name_yrs_attr("ANN")
        expected = " (2000-2002)"

        assert result == expected

    def test_returns_ref_name_and_yrs_averaged_attr_with_climo_dataset_using_short_ref_name(
        self,
    ):
        # Case 1: name is taken from `parameter.short_ref_name`
        param = _create_parameter_object("ref", "climo", self.data_path, "2000", "2002")
        param.short_ref_name = "short_ref_name"
        param.ref_file = self.ref_file

        # Write the climatology dataset out before function call.
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.ref_file}")

        ds1 = Dataset(param, data_type="ref")
        result = ds1.get_name_yrs_attr("ANN")
        expected = "short_ref_name (2000-2002)"

        assert result == expected

    def test_returns_ref_name_and_yrs_averaged_attr_with_climo_dataset_using_reference_name(
        self,
    ):
        # Case 2: name is taken from `parameter.reference_name`
        param = _create_parameter_object("ref", "climo", self.data_path, "2000", "2002")
        param.reference_name = "reference_name"
        param.ref_file = self.ref_file

        # Write the climatology dataset out before function call.
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.ref_file}")

        ds2 = Dataset(param, data_type="ref")
        result = ds2.get_name_yrs_attr("ANN")
        expected = "reference_name (2000-2002)"

        assert result == expected

    def test_returns_ref_name_and_yrs_averaged_attr_with_climo_dataset_using_ref_name(
        self,
    ):
        # Case 3: name is taken from `parameter.ref_name`
        param = _create_parameter_object("ref", "climo", self.data_path, "2000", "2002")
        param.ref_name = "ref_name"
        param.ref_file = self.ref_file

        # Write the climatology dataset out before function call.
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.ref_file}")

        ds3 = Dataset(param, data_type="ref")
        result = ds3.get_name_yrs_attr("ANN")
        expected = "ref_name (2000-2002)"

        assert result == expected

    def test_returns_only_yrs_averaged_attr_if_ref_name_is_not_set_with_climo_dataset(
        self,
    ):
        param = _create_parameter_object("ref", "climo", self.data_path, "2000", "2002")
        param.ref_name = ""
        param.ref_file = self.ref_file

        # Write the climatology dataset out before function call.
        self.ds_climo.to_netcdf(f"{self.data_path}/{param.ref_file}")

        ds3 = Dataset(param, data_type="ref")
        result = ds3.get_name_yrs_attr("ANN")
        expected = " (2000-2002)"

        assert result == expected

    def test_returns_test_name_and_years_averaged_as_single_string_with_timeseries_dataset(
        self,
    ):
        param1 = _create_parameter_object(
            "test", "time_series", self.data_path, "1800", "1850"
        )
        param1.short_test_name = "short_test_name"

        ds1 = Dataset(param1, data_type="test")
        result = ds1.get_name_yrs_attr("ANN")
        expected = "short_test_name (1800-1850)"

        assert result == expected
        assert result == expected
