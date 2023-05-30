from collections import OrderedDict

import pytest

from e3sm_diags.derivations.acme_new import DERIVED_VARIABLES
from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter


class TestDataSetProperties:
    def test_property_is_timeseries_returns_true_and_is_climo_returns_true_for_ref(
        self,
    ):
        parameter = CoreParameter()
        parameter.ref_timeseries_input = True
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="ref")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_true_and_is_climo_returns_false_for_test(
        self,
    ):
        parameter = CoreParameter()
        parameter.test_timeseries_input = True
        parameter.test_data_path = "/path/to/my/files"
        parameter.test_start_yr = "2000"  # type: ignore
        parameter.test_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="test")

        assert ds.is_time_series
        assert not ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_test(
        self,
    ):
        parameter = CoreParameter()
        parameter.ref_timeseries_input = False
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="ref")

        assert not ds.is_time_series
        assert ds.is_climo

    def test_property_is_timeseries_returns_false_and_is_climo_returns_true_for_ref(
        self,
    ):
        parameter = CoreParameter()
        parameter.test_timeseries_input = False
        parameter.test_data_path = "/path/to/my/files"
        parameter.test_start_yr = "2000"  # type: ignore
        parameter.test_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="test")

        assert not ds.is_time_series
        assert ds.is_climo


class TestInit:
    @pytest.fixture(autouse=True)
    def setup(self):
        self.parameter = CoreParameter()

    def test_sets_attrs_if_type_attr_is_ref(self):
        parameter = CoreParameter()
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="ref")

        assert ds.data_path == parameter.reference_data_path
        assert ds.start_yr == parameter.ref_start_yr  # type: ignore
        assert ds.end_yr == parameter.ref_end_yr  # type: ignore

    def test_sets_attrs_if_type_attr_is_test(self):
        parameter = CoreParameter()
        parameter.test_data_path = "/path/to/my/files"
        parameter.test_start_yr = "2000"  # type: ignore
        parameter.test_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="test")

        assert ds.data_path == parameter.test_data_path
        assert ds.start_yr == parameter.test_start_yr  # type: ignore
        assert ds.end_yr == parameter.test_end_yr  # type: ignore

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
        parameter = CoreParameter()
        parameter.sets[0] = "diurnal_cycle"
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="ref")

        assert ds.is_sub_monthly

        parameter.sets[0] = "arm_diags"
        ds2 = Dataset(parameter, type="ref")

        assert ds2.is_sub_monthly

    def test_sets_derived_vars_map(self):
        parameter = CoreParameter()
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore

        ds = Dataset(parameter, type="ref")

        assert ds.derived_vars_map == DERIVED_VARIABLES

    def test_sets_drived_vars_map_with_existing_entry(self):
        parameter = CoreParameter()
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore
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
        parameter = CoreParameter()
        parameter.reference_data_path = "/path/to/my/files"
        parameter.ref_start_yr = "2000"  # type: ignore
        parameter.ref_end_yr = "2001"  # type: ignore
        parameter.derived_variables = {
            "NEW_DERIVED_VAR": OrderedDict([(("some_var",), lambda some_var: some_var)])
        }

        ds = Dataset(parameter, type="ref")

        # The expected `derived_vars_map` result.
        expected = DERIVED_VARIABLES.copy()
        expected["NEW_DERIVED_VAR"] = parameter.derived_variables["NEW_DERIVED_VAR"]  # type: ignore

        assert ds.derived_vars_map == expected


class TestGetClimoDataset:
    def test_raises_error_if_var_arg_is_not_valid(self):
        assert 0

    def test_raises_error_if_season_arg_is_not_valid(self):
        assert 0

    def test_returns_climo_dataset_using_climo_file(self):
        assert 0

    def test_returns_climo_dataset_using_time_series_file(self):
        assert 0

    def test_raises_error_if_invalid_type_specified(self):
        assert 0


class TestGetTimeSeriesDataset:
    def test_raises_error_if_data_is_not_time_series(self):
        assert 0

    def test_raises_error_if_var_arg_is_not_valid(self):
        assert 0

    def test_returns_time_series_dataset(self):
        assert 0

    def test_returns_time_series_dataset_with_centered_time_if_single_point(self):
        assert 0
