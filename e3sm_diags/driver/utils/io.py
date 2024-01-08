from __future__ import annotations

import errno
import json
import os
from typing import Callable

import xarray as xr

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)


def _save_data_metrics_and_plots(
    parameter: CoreParameter,
    plot_func: Callable,
    var_key: str,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
    metrics_dict: MetricsDict | None,
):
    """Save data (optional), metrics, and plots.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.
    plot_func: Callable
        The plot function for the diagnostic set.
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_ref : xr.Dataset | None
        The optional reference dataset. If the diagnostic is a model-only run,
        then it will be None.
    ds_diff : xr.Dataset | None
        The optional difference dataset. If the diagnostic is a model-only run,
        then it will be None.
    metrics_dict : Metrics | None
        The optional dictionary containing metrics for the variable. Some sets
        such as cosp_histogram only calculate spatial average and do not
        use ``metrics_dict``.
    """
    if parameter.save_netcdf:
        _write_vars_to_netcdf(
            parameter,
            var_key,
            ds_test,
            ds_ref,
            ds_diff,
        )

    output_dir = _get_output_dir(parameter)
    filename = f"{parameter.output_file}.json"
    filepath = os.path.join(output_dir, filename)

    if metrics_dict is not None:
        with open(filepath, "w") as outfile:
            json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in {filepath}")

    # Set the viewer description to the "long_name" attr of the variable.
    parameter.viewer_descr[var_key] = ds_test[var_key].attrs.get(
        "long_name", "No long_name attr in test data"
    )

    # Get the function arguments and pass to the set's plotting function.
    args = [
        parameter,
        ds_test[var_key],
        ds_ref[var_key] if ds_ref is not None else None,
        ds_diff[var_key] if ds_diff is not None else None,
    ]
    if metrics_dict is not None:
        args = args + [metrics_dict]

    plot_func(*args)


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
    ds_ref : xr.Dataset
        The dataset containing the ref variable. If this is a model-only run
        then it will be the same dataset as ``ds_test``.
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
