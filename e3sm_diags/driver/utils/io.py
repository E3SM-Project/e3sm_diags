from __future__ import annotations

import errno
import os

import xarray as xr

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)


def _write_vars_to_netcdf(
    parameter: CoreParameter,
    var_key,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
):
    """Saves the test, reference, and difference variables to netCDF files.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter object used to configure the diagnostic runs for the
        sets. The referenced attributes include `save_netcdf, `current_set`,
        `var_id`, `ref_name`, and `output_file`, `results_dir`, and `case_id`.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset | None
        The optional dataset containing the reference variable.
    ds_diff : Optional[xr.DataArray]
        The optional dataset containing the difference between the test and
        reference variables.

    Notes
    -----
    Replaces `e3sm_diags.driver.utils.general.save_ncfiles()`.
    """
    dir_path = _get_output_dir(parameter)
    filename = f"{parameter.output_file}_output.nc"
    output_file = os.path.join(dir_path, filename)

    ds_output = xr.Dataset()
    ds_output[f"{var_key}_test"] = ds_test[var_key]

    if ds_ref is not None:
        ds_output[f"{var_key}_ref"] = ds_ref[var_key]

    if ds_diff is not None:
        ds_output[f"{var_key}_diff"] = ds_diff[var_key]

    ds_output.to_netcdf(output_file)

    logger.info(f"'{var_key}' variable outputs saved to `{output_file}`.")


def _get_output_dir(parameter: CoreParameter):
    """Get the absolute dir path to store the outputs for a diagnostic run.

    If the directory does not exist, attempt to create it.

    When running e3sm_diags is executed with parallelism, a process for another
    set can create the dir already so we can ignore creating the dir for this
    set.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter object used to configure the diagnostic runs for the sets.
        The referenced attributes include `current_set`, `results_dir`, and
        `case_id`.

    Raises
    ------
    OSError
        If the directory does not exist and could not be created.
    """
    results_dir = parameter.results_dir
    dir_path = os.path.join(results_dir, parameter.current_set, parameter.case_id)

    if not os.path.exists(dir_path):
        try:
            os.makedirs(dir_path, 0o755)
        except OSError as e:
            # For parallel runs, raise errors for all cases except when a
            # process already created the directory.
            if e.errno != errno.EEXIST:
                raise OSError(e)

    return dir_path
