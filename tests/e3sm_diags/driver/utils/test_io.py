from copy import deepcopy
from pathlib import Path

import pytest
import xarray as xr

from e3sm_diags.driver.utils.io import _write_vars_to_netcdf
from e3sm_diags.parameter.core_parameter import CoreParameter


class TestWriteVarsToNetcdf:
    @pytest.fixture(autouse=True)
    def setup(self, tmp_path: Path):
        self.param = CoreParameter()
        self.param.save_netcdf = True
        self.param.ref_name = "ref_name"

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
        self.test = xr.DataArray(name="ts", data=[1, 1, 1])
        self.ref = xr.DataArray(name="ts", data=[1, 1, 1])
        self.diff = self.test - self.ref

    def test_does_not_write_to_file_if_parameter_save_netcdf_attr_is_false(self):
        param = CoreParameter()
        param.save_netcdf = False

        _write_vars_to_netcdf(param, self.test, self.ref, diff=None)

        with pytest.raises(FileNotFoundError):
            xr.open_dataarray(f"{self.dir}/ts_test.nc")
            xr.open_dataarray(f"{self.dir}/ts_ref.nc")

    def test_writes_test_variable_to_file(self):
        _write_vars_to_netcdf(self.param, self.test, self.ref, diff=None)

        da_test = xr.open_dataarray(f"{self.dir}/ts_test.nc")
        xr.testing.assert_identical(da_test, self.test)

    def test_writes_ref_variable_to_file_if_param_ref_name_attr_is_set(self):
        _write_vars_to_netcdf(self.param, self.test, self.ref, diff=None)

        da_ref = xr.open_dataarray(f"{self.dir}/ts_ref.nc")
        xr.testing.assert_identical(da_ref, self.ref)

    def test_does_not_write_ref_variable_to_file_if_param_ref_name_attr_is_not_set(
        self,
    ):
        param = deepcopy(self.param)
        param.ref_name = ""

        _write_vars_to_netcdf(param, self.test, self.ref, diff=None)

        with pytest.raises(FileNotFoundError):
            xr.open_dataarray(f"{self.dir}/ts_ref.nc")

    def test_sets_dataarray_names_to_parameter_var_id_attr_if_None(self):
        param = deepcopy(self.param)
        param.var_id = "test_var_id"

        test = self.test.copy()
        ref = self.test.copy()
        diff = self.test.copy()

        test.name = None
        ref.name = None
        diff.name = None

        _write_vars_to_netcdf(param, test, ref, diff=diff)

        da_test = xr.open_dataarray(f"{self.dir}/ts_test.nc")
        da_ref = xr.open_dataarray(f"{self.dir}/ts_ref.nc")
        da_diff = xr.open_dataarray(f"{self.dir}/ts_diff.nc")

        assert da_test.name == param.var_id
        assert da_ref.name == param.var_id
        assert da_diff.name == f"{param.var_id}_diff"
