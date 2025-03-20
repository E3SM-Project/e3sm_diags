from __future__ import annotations

import errno
import json
import os
from typing import Any, Callable, Dict, Literal

import xarray as xr

from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = _setup_child_logger(__name__)


def _save_data_metrics_and_plots(
    parameter: CoreParameter,
    plot_func: Callable,
    var_key: str,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
    metrics_dict: MetricsDict | None,
    plot_kwargs: Dict[str, Any] | None = None,
    viewer_descr: str | None = None,
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
    plot_kwargs : Dict[str, Any] | None
        An optional dictionary containing extra keyword arguments used by a
        plotter, by default None. For example, the enso_diags plotter has extra
        kwargs for confidence levels called `da_test_conf_lvls` and
        `da_ref_conf_lvls`.
    viewer_descr : str | None
        An optional viewer description, by default None. For example,
        the enso_diags driver has a custom viewer description that is not
        the "long_name" variable attribute.
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

    # Set the viewer description to the "long_name" attr of the variable if not
    # manually set.
    if viewer_descr is not None:
        parameter.viewer_descr[var_key] = viewer_descr
    else:
        parameter.viewer_descr[var_key] = ds_test[var_key].attrs.get(
            "long_name", var_key
        )

    # Get the function arguments and pass to the set's plotting function.
    args = {
        "parameter": parameter,
        "da_test": ds_test[var_key],
        "da_ref": ds_ref[var_key] if ds_ref is not None else None,
        "da_diff": ds_diff[var_key] if ds_diff is not None else None,
    }
    if metrics_dict is not None:
        args["metrics_dict"] = metrics_dict

    if plot_kwargs is not None:
        args = {**args, **plot_kwargs}

    plot_func(**args)


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
    _write_to_netcdf(parameter, ds_test[var_key], var_key, "test")

    if ds_ref is not None:
        _write_to_netcdf(parameter, ds_ref[var_key], var_key, "ref")

    if ds_diff is not None:
        _write_to_netcdf(parameter, ds_diff[var_key], var_key, "diff")


def _write_to_netcdf(
    parameter: CoreParameter,
    var: xr.DataArray,
    var_key: str,
    data_type: Literal["test", "ref", "diff"],
):
    filename, filepath = _get_output_filename_filepath(parameter, data_type)

    var.to_netcdf(filepath)

    logger.info(f"'{var_key}' {data_type} variable output saved in: {filepath}")

    return filename


def _get_output_filename_filepath(parameter: CoreParameter, data_type: str):
    dir_path = _get_output_dir(parameter)
    filename = f"{parameter.output_file}_{data_type}.nc"

    filepath = os.path.join(dir_path, filename)

    return filename, filepath


def _write_vars_to_single_netcdf(
    parameter: CoreParameter,
    var_key,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
):
    """Saves the test, reference, and difference variables to a single netCDF.

    NOTE: This function is not currently being used because we need to save
    individual netCDF files (`_write_vars_to_netcdf()`) to perform regression
    testing against the `main` branch, which saves files individually.

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
                raise OSError(e) from e

    return dir_path
