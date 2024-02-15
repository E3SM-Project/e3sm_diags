from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "zonal_mean_xy"
SET_DIR = "debug-654-zonal_mean_xy"
# CFG_PATH = "auxiliary_tools/cdat_regression_testing/654-zonal_mean_xy/debug_zonal_mean_xy_model_vs_obs.cfg"

# run_set(SET_NAME, SET_DIR, CFG_PATH)
run_set(SET_NAME, SET_DIR)


"""
2024-02-15 09:47:22,274 [INFO]: zonal_mean_xy_driver.py(run_diag:68) >> Variable: FLUT
2024-02-15 09:47:22,644 [ERROR]: core_parameter.py(_run_diag:341) >> Error in e3sm_diags.driver.zonal_mean_xy_driver
Traceback (most recent call last):
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py", line 338, in _run_diag
    single_result = module.run_diag(self)
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/zonal_mean_xy_driver.py", line 72, in run_diag
    parameter._set_name_yrs_attrs(test_ds, ref_ds, season)
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py", line 293, in _set_name_yrs_attrs
    self.ref_name_yrs = ds_ref.get_name_yrs_attr(season)
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 166, in get_name_yrs_attr
    diag_name = self._get_ref_name()
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 228, in _get_ref_name
    raise AttributeError(
AttributeError: Either `parameter.short_ref_name`, `parameter.reference_name`, or `parameter.ref_name` must be set to get the name and years attribute for reference datasets.
"""
