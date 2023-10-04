import logging
import os
from copy import deepcopy
from pathlib import Path

import pytest
import xarray as xr

from e3sm_diags.driver.utils.io import _get_output_dir, _write_vars_to_netcdf
from e3sm_diags.parameter.core_parameter import CoreParameter


class TestWriteVarsToNetcdf:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path: Path):
        self.param = CoreParameter()
        self.var_key = "ts"

        # Need to prepend with tmp_path because we use pytest to create temp
        # dirs for storing files temporarily for the test runs.
        self.param.results_dir = f"{tmp_path}/results_dir"
        self.param.current_set = "lat_lon"
        self.param.case_id = "lat_lon_MERRA"
        self.param.output_file = "ts"

        # Create the results directory, which uses the CoreParameter attributes.
        # Example: "<results_dir>/<current_set>/<case_id>/<output_file>_test.nc>"
        self.dir = (
            tmp_path / "results_dir" / self.param.current_set / self.param.case_id
        )
        self.dir.mkdir(parents=True)

        # Input variables for the function
        self.var_key = "ts"
        self.ds_test = xr.Dataset(
            data_vars={"ts": xr.DataArray(name="ts", data=[1, 1, 1])}
        )
        self.ds_ref = xr.Dataset(
            data_vars={"ts": xr.DataArray(name="ts", data=[2, 2, 2])}
        )
        self.ds_diff = self.ds_test - self.ds_ref

    def test_writes_test_variable_to_file(self, caplog):
        # Silence info logger message about saving to a directory.
        caplog.set_level(logging.CRITICAL)

        _write_vars_to_netcdf(self.param, self.var_key, self.ds_test, None, None)

        expected = self.ds_test.copy()
        expected = expected.rename_vars({"ts": "ts_test"})

        result = xr.open_dataset(f"{self.dir}/{self.var_key}_output.nc")
        xr.testing.assert_identical(expected, result)

    def test_writes_ref_and_diff_variables_to_file(self, caplog):
        # Silence info logger message about saving to a directory.
        caplog.set_level(logging.CRITICAL)

        _write_vars_to_netcdf(
            self.param, self.var_key, self.ds_test, self.ds_ref, self.ds_diff
        )

        expected = self.ds_test.copy()
        expected = expected.rename_vars({"ts": "ts_test"})
        expected["ts_ref"] = self.ds_ref["ts"].copy()
        expected["ts_diff"] = self.ds_diff["ts"].copy()

        result = xr.open_dataset(f"{self.dir}/{self.var_key}_output.nc")
        xr.testing.assert_identical(expected, result)


class TestGetOutputDir:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path):
        self.data_path = tmp_path / "input_data"
        self.data_path.mkdir()

        self.param = CoreParameter()
        self.param.results_dir = self.data_path
        self.param.current_set = "lat_lon"
        self.param.case_id = "lat_lon_MERRA"

    def test_raises_error_if_the_directory_does_not_exist_and_cannot_be_created_due_to_permissions(
        self, tmp_path
    ):
        data_path_restricted = tmp_path / "input_data"
        os.chmod(data_path_restricted, 0o444)

        param = deepcopy(self.param)
        param.results_dir = data_path_restricted

        with pytest.raises(OSError):
            _get_output_dir(param)

    def test_creates_directory_if_it_does_not_exist_and_returns_dir_path(self):
        param = CoreParameter()
        param.results_dir = self.data_path
        param.current_set = "lat_lon"
        param.case_id = "lat_lon_MERRA"

        result = _get_output_dir(param)
        assert result == f"{param.results_dir}/{param.current_set}/{param.case_id}"

    def test_ignores_creating_directory_if_it_exists_returns_dir_path(self):
        dir_path = (
            f"{self.param.results_dir}/{self.param.current_set}/{self.param.case_id}"
        )

        nested_dir_path = self.data_path / dir_path
        nested_dir_path.mkdir(parents=True, exist_ok=True)

        result = _get_output_dir(self.param)

        assert result == dir_path
