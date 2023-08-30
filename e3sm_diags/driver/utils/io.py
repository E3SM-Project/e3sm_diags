import errno
import os
from typing import Optional

import xarray as xr

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)


def _write_vars_to_netcdf(
    parameter: CoreParameter,
    test: xr.DataArray,
    ref: xr.DataArray,
    diff: Optional[xr.DataArray],
):
    """Saves the test, reference, and difference variables to netCDF files.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter object used to configure the diagnostic runs for the
        sets. The referenced attributes include `save_netcdf, `current_set`,
        `var_id`, `ref_name`, and `output_file`, `results_dir`, and `case_id`.
    test : xr.DataArray
        The test variable.
    ref : xr.DataArray
        The reference variable.
    diff : Optional[xr.DataArray]
        The optional difference between the test and reference variables.

    Notes
    -----
    This function is intended to replace
    `e3sm_diags.driver.utils.general.save_ncfiles()`.

    # TODO: Do we need to consider the append option like in save_ncfiles()?
    """
    dir_path = _get_output_dir(parameter)
    file_prefix = parameter.output_file
    test_filepath = os.path.join(dir_path, f"{file_prefix}_test.nc")
    ref_filepath = os.path.join(dir_path, f"{file_prefix}_ref.nc")
    diff_filepath = os.path.join(dir_path, f"{file_prefix}_diff.nc")

    if test.name is None:
        test.name = parameter.var_id

    test.to_netcdf(test_filepath)
    logger.info(f"'{test.name}' test variable was saved to: {test_filepath}")

    # Only write out the reference variable if the reference name is set.
    if parameter.ref_name != "":
        if ref.name is None:
            ref.name = parameter.var_id

        ref.to_netcdf(ref_filepath)
        logger.info(f"'{ref.name}' test variable was saved to: {ref_filepath}")

    if diff is not None:
        if diff.name is None:
            diff.name = f"{parameter.var_id}_diff"

        diff.to_netcdf(diff_filepath)
        logger.info(f"'{diff.name}' test variable was saved to: `{diff_filepath}`")


def _get_output_dir(parameter: CoreParameter):
    """Get the absolute directory path to store the outputs for a diagnostic run.

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

    # Create the directory if it does not exist. When running diags in parallel,
    # a process for another set can create the dir already so we can ignore
    # creating the dir for this set.
    if not os.path.exists(dir_path):
        # TODO: Write a unit test for this try/except statement.
        try:
            os.makedirs(dir_path, 0o755)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise OSError(e)

    return dir_path
