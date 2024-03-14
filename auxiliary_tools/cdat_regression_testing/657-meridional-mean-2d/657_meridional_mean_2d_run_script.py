# python -m auxiliary_tools.cdat_regression_testing.792-lat-lon-run-script.792_lat_lon_run_script
from auxiliary_tools.cdat_regression_testing.base_run_script import run_set

SET_NAME = "meridional_mean_2d"
SET_DIR = "657-meridional-mean-2d"
# CFG_PATH: str | None = None
CFG_PATH: str | None = "/global/u2/v/vo13/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/657-meridional-mean-2d/657_meridional_mean_2d.cfg"
MULTIPROCESSING = False

# %%
run_set(SET_NAME, SET_DIR, CFG_PATH, MULTIPROCESSING)


"""
Traceback (most recent call last):
  File "/global/u2/v/vo13/mambaforge/envs/e3sm_diags_dev_nompi_659/lib/python3.10/site-packages/xcdat/regridder/regrid2.py", line 134, in _output_axis_sizes
    axis_name = axis_name_map[standard_name]
KeyError: 'lev'

During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "<string>", line 1, in <module>
  File "/global/u2/v/vo13/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/regrid.py", line 406, in align_grids_to_lower_res
    ds_b_regrid = ds_b_new.regridder.horizontal(
  File "/global/u2/v/vo13/mambaforge/envs/e3sm_diags_dev_nompi_659/lib/python3.10/site-packages/xcdat/regridder/accessor.py", line 324, in horizontal
    output_ds = regridder.horizontal(data_var, self._ds)
  File "/global/u2/v/vo13/mambaforge/envs/e3sm_diags_dev_nompi_659/lib/python3.10/site-packages/xcdat/regridder/regrid2.py", line 96, in horizontal
    output_axis_sizes = self._output_axis_sizes(input_data_var)
  File "/global/u2/v/vo13/mambaforge/envs/e3sm_diags_dev_nompi_659/lib/python3.10/site-packages/xcdat/regridder/regrid2.py", line 136, in _output_axis_sizes
    raise RuntimeError(
RuntimeError: Could not find axis 'lev', ensure 'lev' exists and the attributes are correct.
"""
