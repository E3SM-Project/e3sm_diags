from __future__ import annotations

import collections
import json
import os
from typing import TYPE_CHECKING, Callable, Dict, List, Tuple

import numpy as np
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.diurnal_cycle_xr import composite_diurnal_cycle
from e3sm_diags.driver.utils.io import _get_output_dir
from e3sm_diags.driver.utils.regrid import has_z_axis, regrid_z_axis_to_plevs
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot import arm_diags_plot


if TYPE_CHECKING:
    from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter

logger = _setup_child_logger(__name__)

RefsTestMetrics = collections.namedtuple(
    "RefsTestMetrics", ["refs", "test", "metrics", "misc"]
)

# A dictionary that maps diags_set to the appropriate seasons for grouping.
SEASONS_BY_DIAG: Dict[str, List[ClimoFreq]] = {
    "diurnal_cycle": ["DJF", "MAM", "JJA", "SON"],
    "annual_cycle": ["ANNUALCYCLE"],
    "diurnal_cycle_zt": ["ANNUALCYCLE"],
}


def run_diag(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    """Run the specified diagnostic set based on the given ARMDiagsParameter.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The ARMDiagsParameter object containing the configuration for the diagnostic run.

    Returns
    -------
    ARMDiagsParameter
        The updated ARMDiagsParameter object after running the diagnostic set.

    Raises
    ------
    RuntimeError
        If the specified diags_set is invalid.
    """
    if parameter.diags_set == "annual_cycle":
        return _run_diag_annual_cycle(parameter)
    elif parameter.diags_set == "diurnal_cycle":
        return _run_diag_diurnal_cycle(parameter)
    elif parameter.diags_set == "diurnal_cycle_zt":
        return _run_diag_diurnal_cycle_zt(parameter)
    elif parameter.diags_set == "pdf_daily":
        logger.info("'run_diag_pdf_daily' is not yet implemented.")
    elif parameter.diags_set == "convection_onset":
        return _run_diag_convection_onset(parameter)
    elif parameter.diags_set == "aerosol_activation":
        return _run_diag_aerosol_activation(parameter)

    raise RuntimeError(f"Invalid diags_set={parameter.diags_set}")


def _run_diag_diurnal_cycle(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = SEASONS_BY_DIAG["diurnal_cycle"]
    test_ds = Dataset(parameter, data_type="test")

    for region in regions:
        logger.info(f"Selected region: {region}")
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info(f"Season: {season}")

            for var in variables:
                logger.info(f"Variable: {var}")

                ds_test = test_ds.get_time_series_dataset(var, single_point=True)
                test_diurnal, lst = composite_diurnal_cycle(  # type: ignore
                    ds_test, var, season, fft=False
                )

                refs = []

                if "armdiags" in ref_name:
                    if region != "sgpc1":
                        raise RuntimeError(
                            f"Diurnal cycle of {region} at Site: {var} is not "
                            "supported yet"
                        )
                    else:
                        ref_file_name = "sgparmdiagsmondiurnalC1.c1.nc"
                        ref_file = os.path.join(ref_path, ref_file_name)
                        ds_ref = xr.open_dataset(ref_file)

                        if var == "PRECT":
                            # Converting mm/second to mm/day"
                            ds_ref[var] = ds_ref["pr"] * 3600.0 * 24
                            ds_ref["lat"] = ds_test.lat.values
                            ds_ref["lon"] = ds_test.lon.values

                            ref_diurnal, lst = composite_diurnal_cycle(  # type: ignore
                                ds_ref, var, season, fft=False
                            )

                            ref = ref_diurnal

                else:
                    ref_data = Dataset(parameter, data_type="ref")
                    ds_ref = ref_data.get_time_series_dataset(var, single_point=True)

                    ref_diurnal, lst = composite_diurnal_cycle(  # type: ignore
                        ds_ref, var, season, fft=False
                    )
                    ref = ref_diurnal

                refs.append(ref)

                # Create the metrics dictionary.
                metrics_dict = {}
                metrics_dict["unit"] = ds_test[var].units

                # Update the result metrics dictionary and store it in the vars_
                # to_data dictionary.
                result = RefsTestMetrics(
                    test=test_diurnal, refs=refs, metrics=None, misc=lst
                )
                vars_to_data[season] = result

                parameter.output_file = "-".join([ref_name, var, season, region])
                _save_metrics_to_json(parameter, metrics_dict)

                # Set the plot and viewer output attributes.
                parameter.viewer_descr[var] = ds_test[var].attrs.get("long_name", var)
                parameter.test_name_yrs = test_ds.get_name_yrs_attr()
                parameter.var_name = ds_test[var].attrs.get("long_name", var)
                parameter.var_units = ds_test[var].attrs.get("units", var)

                arm_diags_plot._plot_diurnal_cycle(parameter, vars_to_data[season])

    return parameter


def _run_diag_diurnal_cycle_zt(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = SEASONS_BY_DIAG["diurnal_cycle_zt"]
    plevs: List[float] = list(np.linspace(100, 1000, 37))

    for region in regions:
        logger.info(f"Selected region: {region}")
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info(f"Season: {season}")
            for var in variables:
                logger.info(f"Variable: {var}")

                test_ds = Dataset(parameter, data_type="test")
                ds_test = test_ds.get_time_series_dataset(var, single_point=True)

                if has_z_axis(ds_test[var]):
                    ds_test_plevs = regrid_z_axis_to_plevs(ds_test, var, plevs)
                    ds_test_plevs["lat"] = ds_test.lat
                    ds_test_plevs["lon"] = ds_test.lon

                    test_diurnal, lst = composite_diurnal_cycle(  # type: ignore
                        ds_test_plevs, var, season, fft=False
                    )

                refs = []
                if "armdiags" in ref_name:
                    ref_file_name = (
                        region[:3]
                        + "armdiagsmondiurnalclim"
                        + region[3:5].upper()
                        + ".c1.nc"
                    )

                    ref_file = os.path.join(ref_path, ref_file_name)
                    ds_ref = xr.open_dataset(ref_file)

                    if var == "CLOUD":
                        ref_var = ds_ref["cl_p"].values
                        ref_var = np.reshape(ref_var, (12, 24, ref_var.shape[1]))
                        ref_diurnal = ref_var

                else:
                    ref_ds = Dataset(parameter, data_type="ref")
                    ds_ref = ref_ds.get_time_series_dataset(var, single_point=True)

                    ds_ref_plevs = regrid_z_axis_to_plevs(ds_ref, var, plevs)
                    ds_ref_plevs["lat"] = ds_test.lat.values
                    ds_ref_plevs["lon"] = ds_test.lon.values

                    ref_diurnal, lst = composite_diurnal_cycle(  # type: ignore
                        ds_ref_plevs, var, season, fft=False
                    )

                refs.append(ref_diurnal)

                # Create the metrics dictionary.
                metrics_dict = {}
                metrics_dict["unit"] = ds_test[var].units

                # Update the result metrics dictionary and store it in the vars_
                # to_data dictionary.
                result = RefsTestMetrics(
                    test=test_diurnal, refs=refs, metrics=None, misc=lst
                )
                vars_to_data[season] = result

                # Save the metrics to json.
                parameter.output_file = "-".join([ref_name, var, season, region])
                _save_metrics_to_json(parameter, metrics_dict)

                # Save the plot and viewer output attributes.
                parameter.viewer_descr[var] = ds_test[var].attrs.get("long_name", var)
                parameter.test_name_yrs = test_ds.get_name_yrs_attr()
                parameter.var_name = ds_test[var].attrs.get("long_name", var)
                parameter.var_units = ds_test[var].attrs.get("units", var)

                arm_diags_plot._plot_diurnal_cycle_zt(parameter, vars_to_data[season])

    return parameter


def _run_diag_annual_cycle(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = SEASONS_BY_DIAG["annual_cycle"]

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info(f"Selected region: {region}")
        vars_to_data = collections.OrderedDict()

        test_ds = Dataset(parameter, data_type="test")

        for season in seasons:
            logger.info(f"Season: {season}")

            for var in variables:
                logger.info(f"Variable: {var}")

                ds_test = test_ds.get_climo_dataset(var, season)
                da_test = ds_test[var]

                refs = []

                if "armdiags" in ref_name:
                    ref_file = os.path.join(
                        ref_path,
                        region[:3] + "armdiagsmon" + region[3:5].upper() + ".c1.nc",
                    )

                    ds_ref = xr.open_dataset(ref_file)
                    vars_funcs = _get_vars_funcs_for_derived_var(ds_ref, var)
                    target_var = list(vars_funcs.keys())[0][0]

                    # NOTE: The bounds dimension can be "nv", which is not
                    # currently recognized as a valid bounds dimension
                    # by xcdat. We rename it to "bnds" to make it compatible.
                    ds_ref = _rename_bounds_dim(ds_ref)

                    ds_ref_climo = ds_ref.temporal.climatology(target_var, "month")
                    da_ref = vars_funcs[(target_var,)](ds_ref_climo[target_var]).rename(
                        var
                    )
                    if da_ref.attrs.get("standard_name") is not None:
                        da_ref.attrs["long_name"] = da_ref.attrs["standard_name"]
                else:
                    ref_ds = Dataset(parameter, data_type="ref")
                    ds_ref = ref_ds.get_climo_dataset(var, season)
                    da_ref = ds_ref[var]

                # TODO make this module work with global monthly data
                # ref_domain = utils.regrid._subset_on_arm_coord(ref, var, region)
                ref_domain = da_ref.values
                # ref[var].ref_name = ref_name
                refs.append(ref_domain)

                # TODO make this module work with global monthly data
                # test_domain = utils.regrid._subset_on_arm_coord(test, var, region)
                test_domain = da_test.values

                # Create the metrics dictionary.
                metrics_dict = _get_metrics_dict(test_domain, ref_domain)

                result = RefsTestMetrics(
                    test=test_domain, refs=refs, metrics=metrics_dict, misc=None
                )
                vars_to_data[season] = result
                metrics_dict["unit"] = da_test.units
                metrics_dict["ref_domain"] = list(ref_domain)  # type: ignore
                metrics_dict["test_domain"] = list(test_domain)  # type: ignore

                # Save the metrics to json.
                parameter.output_file = "-".join([ref_name, var, season, region])
                _save_metrics_to_json(parameter, metrics_dict)

                # Set the plot and viewer output attributes.
                parameter.viewer_descr[var] = da_test.attrs.get("long_name", var)
                parameter.test_name_yrs = test_ds.get_name_yrs_attr()
                parameter.var_name = da_test.attrs.get("long_name", var)
                parameter.var_units = da_test.attrs.get("units", var)

            if season == "ANNUALCYCLE":
                arm_diags_plot._plot_annual_cycle(parameter, var, vars_to_data[season])

    return parameter


def _run_diag_convection_onset(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    test_ds = Dataset(parameter, data_type="test")

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info(f"Selected region: {region}")

        ds_test_pr = test_ds.get_time_series_dataset("PRECT", single_point=True)
        test_pr = ds_test_pr["PRECT"].values / 24  # convert to mm/hr

        ds_test_prw = test_ds.get_time_series_dataset("TMQ", single_point=True)
        test_prw = ds_test_prw["TMQ"].values

        if "armdiags" in ref_name:
            if region == "sgp":
                ref_file_name = "sgparmdiags1hrC1.c1.nc"
            else:
                ref_file_name = (
                    region[:3] + "armdiags1hr" + region[3:5].upper() + ".c1.nc"
                )

            ref_file = os.path.join(ref_path, ref_file_name)
            ds_ref = xr.open_dataset(ref_file)

            ref_pr = ds_ref["pr"].values  # mm/hr
            ref_pr[ref_pr < -900] = np.nan

            ref_prw = ds_ref["prw"].values  # mm
            ref_prw[ref_prw < -900] = np.nan
        else:
            ref_ds = Dataset(parameter, data_type="ref")

            ds_ref_pr = ref_ds.get_time_series_dataset("PRECT", single_point=True)
            ref_pr = ds_ref_pr["PRECT"].values / 24

            ds_ref_prw = ref_ds.get_time_series_dataset("TMQ", single_point=True)
            ref_prw = ds_ref_prw["TMQ"].values

        # Set the plot and viewer output attributes.
        parameter.test_name_yrs = test_ds.get_name_yrs_attr()
        parameter.output_file = "-".join([ref_name, "convection-onset", region])

        time_coords = xc.get_dim_coords(ds_test_pr, axis="T")
        parameter.time_interval = int(time_coords[1].dt.hour - time_coords[0].dt.hour)

        arm_diags_plot._plot_convection_onset_statistics(
            parameter, region, test_pr, test_prw, ref_pr, ref_prw
        )

    return parameter


def _run_diag_aerosol_activation(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    # Supported regions are in `e3sm_diags/derivations/default_regions_xr.py`
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    # Possible variables are ccn01, ccn02, ccn05
    variables = parameter.variables

    for region in regions:
        logger.info(f"Selected region: {region}")

        for variable in variables:
            test_data = Dataset(parameter, data_type="test")

            ds_test_a_num = test_data.get_time_series_dataset(
                "a_num", single_point=True
            )
            test_a_num = ds_test_a_num["a_num"].values[:, -1]

            ds_test_ccn = test_data.get_time_series_dataset(variable, single_point=True)
            test_ccn = ds_test_ccn[variable].values[:, -1]

            # Get the name of the data, appended with the years averaged.
            parameter.test_name_yrs = test_data.get_name_yrs_attr()

            if "armdiags" in ref_name:
                ref_file = os.path.join(
                    ref_path,
                    region[:3] + "armdiagsaciactivate" + region[3:5].upper() + ".c1.nc",
                )
                ds_ref = xr.open_dataset(ref_file)

                ref_a_num = ds_ref["cpc_bulk"].values
                ref_ccn = ds_ref[f"{variable}_bulk"].values

            else:
                ref_ds = Dataset(parameter, data_type="test")

                ds_ref_a_num = ref_ds.get_time_series_dataset(
                    "a_num", single_point=True
                )
                ref_a_num = ds_ref_a_num["a_num"].values[:, -1]

                ds_ref_ccn = ref_ds.get_time_series_dataset(variable, single_point=True)
                ref_ccn = ds_ref_ccn[variable].values[:, -1]

            parameter.output_file = "-".join(
                [ref_name, "aerosol-activation", region, variable]
            )
            arm_diags_plot._plot_aerosol_activation(
                parameter, region, variable, test_a_num, test_ccn, ref_a_num, ref_ccn
            )

    return parameter


def _get_vars_funcs_for_derived_var(
    ds: xr.Dataset, var: str
) -> Dict[Tuple[str], Callable]:
    """
    Get a dictionary that maps the list of variables to the derivation function.

    The ARM Diags reference datasets file names do not follow E3SM naming
    convention, this function is a simplified derived variable routine for
    accomodating files from ARM Diags.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    var : str
        The key of the variable.

    Returns
    -------
    Dict[Tuple[str], Callable]
        A tuple of the derived variable and the function to calculate it.
    """
    vars_to_func_dict = DERIVED_VARIABLES[var]
    vars_in_file = set(ds.keys())

    # e.g,. [('pr',), ('PRECC', 'PRECL')]
    possible_vars: List[Tuple[str]] = list(vars_to_func_dict.keys())  # type: ignore

    for list_of_vars in possible_vars:
        if vars_in_file.issuperset(list_of_vars):
            return {list_of_vars: vars_to_func_dict[list_of_vars]}

    raise RuntimeError(f"No derived variable function found for the variable {var}")


def _get_metrics_dict(
    test_var: xr.DataArray | np.ndarray, ref_var: xr.DataArray | np.ndarray
) -> Dict[str, float]:
    """Calculate various metrics between the test and reference DataArrays.

    Parameters
    ----------
    test_var : xr.DataArray | np.ndarray
        The test variable.
    ref_var : xr.DataArray | np.ndarray
        The reference variable.

    Returns
    -------
    metrics_dict : Dict[str, float]
        A dictionary containing the calculated metrics:
        - 'test_mean': The mean of the test_var.
        - 'ref_mean': The mean of the ref_var.
        - 'test_std': The standard deviation of the test_var.
        - 'ref_std': The standard deviation of the ref_var.
        - 'rmse': The root mean squared error between the test and ref vars.
        - 'corr': The correlation coefficient between the test and ref vars.
    """
    return {
        "test_mean": float(test_var.mean()),
        "ref_mean": float(ref_var.mean()),
        "test_std": float(test_var.std(ddof=1)),
        "ref_std": float(ref_var.std(ddof=1)),
        "rmse": float(_rmse(test_var, ref_var)),
        "corr": float(np.corrcoef(test_var, ref_var)[0, 1]),
    }


def _rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def _save_metrics_to_json(parameter: ARMDiagsParameter, metrics_dict: Dict[str, float]):
    """Save metrics dictionary to a JSON file.

    Parameters
    ----------
    parameter : ARMDiagsParameter
        The parameter object.
    metrics_dict : dict
        Dictionary containing the metrics to be saved.
    """
    output_dir = _get_output_dir(parameter)
    filename = parameter.output_file + ".json"
    abs_path = os.path.join(output_dir, filename)

    with open(abs_path, "w") as outfile:
        json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in: {abs_path}")


def _rename_bounds_dim(ds: xr.Dataset) -> xr.Dataset:
    """
    Renames the bounds dimension "nv" to "bnds" in the given xarray.Dataset for
    xCDAT compatibility.

    This is a temporary workaround to ensure compatibility with xCDAT's bounds
    handling. The bounds dimension "nv" is commonly used in datasets to
    represent the number of vertices in a polygon, but xCDAT expects the
    bounds dimension to be in `xcdat.bounds.VALID_BOUNDS_DIMS`. This function
    renames "nv" to "bnds" to align with xCDAT's expectations.

    Parameters
    ----------
    ds : xr.Dataset
        The input xarray.Dataset which may contain a bounds dimension named "nv".

    Returns
    -------
    xr.Dataset
        A new xarray.Dataset with the "nv" dimension renamed to "bnds" if it
        existed; otherwise, the original dataset copy.
    """
    ds_new = ds.copy()

    if "nv" in ds_new.dims:
        ds_new = ds_new.rename({"nv": "bnds"})

    return ds_new
