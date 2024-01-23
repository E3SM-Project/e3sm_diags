"""
The template run script used for generating results on your development branch.

Steps:
1. Activate your conda dev env for your branch
2. `make install` to install the latest version of your branch code into the env
3. Copy this script into `auxiliary_tools/cdat_regression_testing/<ISSUE>-<SET_NAME>`
4. Update `SET_DIR` string variable
5. Update `SET_NAME` string variable.
   - Options include: "lat_lon", "zonal_mean_xy", "zonal_mean_2d",
     "zonal_mean_2d_stratosphere", "polar", "cosp_histogram",
     "meridional_mean_2d", "annual_cycle_zonal_mean", "enso_diags", "qbo",
     "area_mean_time_series", "diurnal_cycle", "streamflow", "arm_diags",
     "tc_analysis", "aerosol_aeronet", "aerosol_budget", "mp_partition",
6. Run this script
   - Make sure to run this command on NERSC perlmutter cpu:
    `salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm
    conda activate <NAME-OF-DEV-ENV>`
   - python auxiliary_tools/cdat_regression_testing/<ISSUE-<SET_NAME>
7. Make a copy of the CDAT regression testing notebook in the same directory
   as this script and follow the instructions there to start testing.
"""
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "cosp_histogram"
SET_DIR = "660-cosp-histogram"

run_set(SET_NAME, SET_DIR)


"""Remaining issues 1-23-24
2024-01-23 15:39:21,094 [INFO]: cosp_histogram_driver.py(run_diag:49) >> Variable: COSP_HISTOGRAM_MODIS
2024-01-23 15:39:46,948 [ERROR]: core_parameter.py(_run_diag:331) >> Error in e3sm_diags.driver.cosp_histogram_driver
Traceback (most recent call last):
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py", line 328, in _run_diag
    single_result = module.run_diag(self)
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/cosp_histogram_driver.py", line 55, in run_diag
    ds_test = test_ds.get_climo_dataset(var_key, season)
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 365, in get_climo_dataset
    ds = self._get_climo_dataset(season)
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 397, in _get_climo_dataset
    ds = self._get_dataset_with_derived_climo_var(ds)
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 638, in _get_dataset_with_derived_climo_var
    ds_final = self._get_dataset_with_derivation_func(
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py", line 831, in _get_dataset_with_derivation_func
    ds_final = func(*func_args)
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/derivations/formulas_cosp.py", line 111, in cosp_histogram_standardize
    prs = _get_cloud_axis(var, "prs")
  File "/global/homes/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/derivations/formulas_cosp.py", line 202, in _get_cloud_axis
    raise KeyError(
KeyError: "The 'prs' axis is not in the 'CLMODIS' to standardize the cosp histogram."
"""
