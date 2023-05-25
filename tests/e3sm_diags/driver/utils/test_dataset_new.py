import pytest

from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.parameter.core_parameter import CoreParameter


class TestDataset:
    @pytest.fixture(autouse=True)
    def setup(self):
        parameter = CoreParameter()
        self.ds = Dataset(parameter, type="ref")

    def test__init__(self):
        parameter = CoreParameter()
        Dataset(parameter, type="ref")

    def test_property_is_timeseries_returns_true(self):
        assert 0

    def test_property_is_timeseries_returns_false(self):
        assert 0

    def test_property_is_climo_returns_true(self):
        assert 0

    def test_property_is_climo_returns_false(self):
        assert 0

    def test_get_static_variable_returns_time_series_variable(self):
        assert 0

    def test_get_static_variable_raises_error_if_dir_has_multiple_time_series_files_for_var(
        self,
    ):
        assert 0

    def test_get_static_variable_raises_error_if_dir_has_multiple_time_series_files_with_ref_name_prepended_for_var(
        self,
    ):
        assert 0

    def test_get_attr_from_climo_returns_attr(self):
        assert 0

    def test_get_attr_from_climo_raises_value_error_if_var_not_found_in_dataset_file(
        self,
    ):
        assert 0

    def test_get_timeseries_variable_returns_timeseries_variable_for_a_given_season_using_derived_vars_dict(
        self,
    ):
        assert 0

    def test_get_timeseries_variable_returns_timeseries_variable_for_a_given_season_if_time_series_file_exists(
        self,
    ):
        assert 0

    def test_get_timeseries_variable_raises_error_if_invalid_year_range_for_time_series_file(
        self,
    ):
        # If start year < var_start_year
        # If end_year > var_start_year
        assert 0

    def test_get_timeseries_variable_raises_error_if_var_file_not_found_or_var_could_not_be_derived(
        self,
    ):
        assert 0

    # get_climo_variable() tests
    # --------------------------
    def test_get_climo_variable_raises_error_if_var_parameter_is_invalid(self):
        assert 0

    def test_get_climo_variable_raises_error_if_season_parameter_is_invalid(self):
        assert 0

    def test_get_climo_variable_returns_ref_climo_var_from_time_series_file(self):
        assert 0

    def test_get_climo_variable_returns_test_climo_var_from_time_series_file(self):
        assert 0

    def test_get_climo_variable_returns_ref_climo_var_from_climo_file(self):
        assert 0

    def test_get_climo_variable_returns_test_climo_var_from_climo_file(self):
        assert 0

    def test_get_climo_variable_raises_error_if_if_unable_to_determine_type_and_source_of_var(
        self,
    ):
        # For example, whether it is a ref or test variable and whether it comes
        # from climatology or time series files.
        assert 0
