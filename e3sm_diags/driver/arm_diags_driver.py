from __future__ import annotations

import collections
import json
import os
from typing import TYPE_CHECKING

import cdms2
import numpy as np
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES
from e3sm_diags.driver import utils
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.diurnal_cycle_xr import composite_diurnal_cycle
from e3sm_diags.driver.utils.regrid import regrid_z_axis_to_plevs
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot import arm_diags_plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter

logger = custom_logger(__name__)

RefsTestMetrics = collections.namedtuple(
    "RefsTestMetrics", ["refs", "test", "metrics", "misc"]
)


def get_vars_funcs_for_derived_var(data_file, var):
    """
    The ARM Diags reference datasets file names do not follow E3SM naming convention, this function is a simplified derived variable routine for accomodating files from ARM Diags.
    """
    vars_to_func_dict = DERIVED_VARIABLES[var]
    vars_in_file = set(data_file.keys())
    # ex: [('pr',), ('PRECC', 'PRECL')]
    possible_vars = list(vars_to_func_dict.keys())  # type: ignore

    for list_of_vars in possible_vars:
        if vars_in_file.issuperset(list_of_vars):
            # All of the variables (list_of_vars) are in data_file.
            # Return the corresponding dict.
            return {list_of_vars: vars_to_func_dict[list_of_vars]}  # type: ignore


def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())


def create_metrics(test, ref):
    """
    For this plotset, calculate the mean, std of test (array-like) and ref (array-like), as well as rmse and corr of two datasets and return a dict of that.
    """
    return {
        "test_mean": float(test.mean()),
        "ref_mean": float(ref.mean()),
        "test_std": float(test.std(ddof=1)),
        "ref_std": float(ref.std(ddof=1)),
        "rmse": float(rmse(test, ref)),
        "corr": float(np.corrcoef(test, ref)[0, 1]),
    }


def run_diag_diurnal_cycle(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = ["DJF", "MAM", "JJA", "SON"]

    for region in regions:
        logger.info("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info("Season: {}".format(season))
            for var in variables:
                logger.info("Variable: {}".format(var))

                test_data = Dataset(parameter, data_type="test")
                # test is a dataset
                test = test_data.get_time_series_dataset(var, single_point=True)
                test_diurnal, lst = composite_diurnal_cycle(
                    test, var, season, fft=False
                )

                parameter.viewer_descr[var] = test[var].long_name
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = test_data.get_name_yrs_attr()
                parameter.var_name = test[var].long_name
                parameter.var_units = test[var].units

                refs = []

                if "armdiags" in ref_name:
                    if region != "sgpc1":
                        msg = "Diurnal cycle of {} at Site: {} is not supported yet".format(
                            region, var
                        )
                        raise RuntimeError(msg)
                    else:
                        ref_file_name = "sgparmdiagsmondiurnalC1.c1.nc"

                        ref_file = os.path.join(ref_path, ref_file_name)
                        ref = xr.open_dataset(ref_file)
                        if var == "PRECT":
                            ref[var] = (
                                ref["pr"] * 3600.0 * 24
                            )  # Converting mm/second to mm/day"
                            ref["lat"] = test.lat.values
                            ref["lon"] = test.lon.values
                            ref_diurnal, lst = composite_diurnal_cycle(
                                ref, var, season, fft=False
                            )

                            ref = ref_diurnal

                else:
                    ref_data = Dataset(parameter, data_type="ref")
                    # ref is a dataset
                    ref = ref_data.get_time_series_dataset(var, single_point=True)
                    ref_diurnal, lst = composite_diurnal_cycle(
                        ref, var, season, fft=False
                    )

                refs.append(ref)

                metrics_dict = {}
                result = RefsTestMetrics(
                    test=test_diurnal, refs=refs, metrics=None, misc=lst
                )
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict["unit"] = test[var].units
                parameter.output_file = "-".join([ref_name, var, season, region])
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                with open(fnm, "w") as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                logger.info("Metrics saved in: " + fnm)

                arm_diags_plot.plot_diurnal_cycle(var, vars_to_data[season], parameter)

    return parameter


def run_diag_diurnal_cycle_zt(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = ["ANNUALCYCLE"]
    plevs = np.linspace(100, 1000, 37)

    for region in regions:
        logger.info("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info("Season: {}".format(season))
            for var in variables:
                logger.info("Variable: {}".format(var))

                test_data = Dataset(parameter, data_type="test")
                # test is a dataset
                test = test_data.get_time_series_dataset(var, single_point=True)
                test_p = regrid_z_axis_to_plevs(test, var, plevs)
                test_p["lat"] = test.lat.values
                test_p["lon"] = test.lon.values
                test_diurnal, lst = composite_diurnal_cycle(
                    test_p, var, season, fft=False
                )

                parameter.viewer_descr[var] = test[var].long_name
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = test_data.get_name_yrs_attr()
                parameter.var_name = test[var].long_name
                parameter.var_units = test[var].units

                refs = []

                if "armdiags" in ref_name:
                    ref_file_name = (
                        region[:3]
                        + "armdiagsmondiurnalclim"
                        + region[3:5].upper()
                        + ".c1.nc"
                    )

                    ref_file = os.path.join(ref_path, ref_file_name)
                    ref = xr.open_dataset(ref_file)
                    if var == "CLOUD":
                        ref_var = ref["cl_p"].values
                        ref_var = np.reshape(ref_var, (12, 24, ref_var.shape[1]))
                        ref_diurnal = ref_var

                else:
                    ref_data = Dataset(parameter, data_type="ref")
                    # test is a dataset
                    ref = ref_data.get_time_series_dataset(var, single_point=True)
                    ref_p = regrid_z_axis_to_plevs(ref, var, plevs)
                    ref_p["lat"] = test.lat.values
                    ref_p["lon"] = test.lon.values
                    ref_diurnal, lst = composite_diurnal_cycle(
                        ref_p, var, season, fft=False
                    )

                refs.append(ref_diurnal)

                metrics_dict = {}
                result = RefsTestMetrics(
                    test=test_diurnal, refs=refs, metrics=None, misc=lst
                )
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict["unit"] = test[var].units
                parameter.output_file = "-".join([ref_name, var, season, region])
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                with open(fnm, "w") as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                logger.info("Metrics saved in: " + fnm)

                arm_diags_plot.plot_diurnal_cycle_zt(
                    var, vars_to_data[season], parameter
                )

    return parameter


def run_diag_annual_cycle(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = ["ANNUALCYCLE"]

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info("Season: {}".format(season))
            for var in variables:
                logger.info("Variable: {}".format(var))
                # CDAT code
                # test_data = utils.dataset.Dataset(parameter, test=True)
                # test = test_data.get_climo_variable(var, season)
                test_data = Dataset(parameter, data_type="test")
                # test is a dataarray
                test = test_data.get_climo_dataset(var, season)[var]

                parameter.viewer_descr[var] = test.long_name
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = test_data.get_name_yrs_attr()
                parameter.var_name = test.long_name
                parameter.var_units = test.units

                refs = []

                if "armdiags" in ref_name:
                    ref_file = os.path.join(
                        ref_path,
                        region[:3] + "armdiagsmon" + region[3:5].upper() + ".c1.nc",
                    )

                    ref_data = xr.open_dataset(ref_file)
                    vars_funcs = get_vars_funcs_for_derived_var(ref_data, var)
                    target_var = list(vars_funcs.keys())[0][0]

                    ref = ref_data.temporal.climatology(target_var, "month")
                    # ref is a dataarray
                    ref = vars_funcs[(target_var,)](ref[target_var]).rename(var)

                    if hasattr(ref, "standard_name"):
                        ref.long_name = ref.standard_name
                else:
                    ref_data = Dataset(parameter, data_type="ref")
                    ref = ref_data.get_climo_dataset(var, season)[var]
                # TODO make this module work with global monthly data
                # ref_domain = utils.regrid._subset_on_arm_coord(ref, var, region)
                ref_domain = ref.values
                # ref[var].ref_name = ref_name
                refs.append(ref_domain)

                # TODO make this module work with global monthly data
                # test_domain = utils.regrid._subset_on_arm_coord(test, var, region)
                test_domain = test.values

                metrics_dict = create_metrics(test_domain, ref_domain)

                result = RefsTestMetrics(
                    test=test_domain, refs=refs, metrics=metrics_dict, misc=None
                )
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict["unit"] = test.units
                metrics_dict["ref_domain"] = list(ref_domain)
                metrics_dict["test_domain"] = list(test_domain)
                parameter.output_file = "-".join([ref_name, var, season, region])
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                with open(fnm, "w") as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                logger.info(f"Metrics saved in: {fnm}")

            if season == "ANNUALCYCLE":
                arm_diags_plot.plot_annual_cycle(var, vars_to_data[season], parameter)

    return parameter


def run_diag_convection_onset(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path
    # Read in observation data

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info("Selected region: {}".format(region))

        test_data = Dataset(parameter, data_type="test")
        test_pr = test_data.get_time_series_dataset("PRECT", single_point=True)
        test_time_coord = test_pr.time
        test_pr = test_pr["PRECT"].values / 24  # convert to mm/hr
        test_prw = test_data.get_time_series_dataset("TMQ", single_point=True)
        test_prw = test_prw["TMQ"].values

        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = test_data.get_name_yrs_attr()

        if "armdiags" in ref_name:
            if region == "sgp":
                ref_file_name = "sgparmdiags1hrC1.c1.nc"
            else:
                ref_file_name = (
                    region[:3] + "armdiags1hr" + region[3:5].upper() + ".c1.nc"
                )
            ref_file = os.path.join(ref_path, ref_file_name)
            ref_data = xr.open_dataset(ref_file)
            # ref_data = cdms2.open(ref_file)
            ref_pr = ref_data["pr"].values  # mm/hr
            ref_pr[ref_pr < -900] = np.nan
            ref_prw = ref_data["prw"].values  # mm
            ref_prw[ref_prw < -900] = np.nan
        else:
            ref_data = Dataset(parameter, data_type="ref")
            ref_pr = test_data.get_time_series_dataset("PRECT", single_point=True)
            ref_pr = ref_pr["PRECT"].values / 24
            ref_prw = test_data.get_time_series_dataset("TMQ", single_point=True)
            ref_prw = ref_prw["TMQ"].values
        parameter.output_file = "-".join([ref_name, "convection-onset", region])
        parameter.time_interval = int(
            test_time_coord[1].dt.hour - test_time_coord[0].dt.hour
        )
        arm_diags_plot.plot_convection_onset_statistics(
            test_pr, test_prw, ref_pr, ref_prw, parameter, region
        )

    return parameter


def run_diag_aerosol_activation(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path
    variables = parameter.variables
    # Read in observation data

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info("Selected region: {}".format(region))
        # Possible variables are ccn01, ccn02, ccn05
        for variable in variables:
            test_data = Dataset(parameter, data_type="test")
            test_a_num = test_data.get_time_series_dataset("a_num", single_point=True)
            test_ccn = test_data.get_time_series_dataset(variable, single_point=True)
            test_a_num = test_a_num["a_num"].values[:, -1]
            test_ccn = test_ccn[variable].values[:, -1]

            # Get the name of the data, appended with the years averaged.
            parameter.test_name_yrs = test_data.get_name_yrs_attr()

            if "armdiags" in ref_name:
                ref_file = os.path.join(
                    ref_path,
                    region[:3] + "armdiagsaciactivate" + region[3:5].upper() + ".c1.nc",
                )
                ref_data = xr.open_dataset(ref_file)
                ref_a_num = ref_data["cpc_bulk"].values
                ref_ccn = ref_data[f"{variable}_bulk"].values

            else:
                ref_data = Dataset(parameter, data_type="test")
                ref_a_num = ref_data.get_time_series_dataset("a_num", single_point=True)
                ref_ccn = ref_data.get_time_series_dataset(variable, single_point=True)
                ref_a_num = ref_a_num["a_num"].values[:, -1]
                ref_ccn = ref_ccn[variable].values[:, -1]

            parameter.output_file = "-".join(
                [ref_name, "aerosol-activation", region, variable]
            )
            arm_diags_plot.plot_aerosol_activation(
                test_a_num, test_ccn, ref_a_num, ref_ccn, parameter, region, variable
            )

    return parameter


def run_diag_annual_cycle_aerosol(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    variables = parameter.variables
    regions = parameter.regions
    ref_name = parameter.ref_name
    ref_path = parameter.reference_data_path

    seasons = ["ANNUALCYCLE"]

    for region in regions:
        # The regions that are supported are in e3sm_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        logger.info("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            logger.info("Season: {}".format(season))
            for var in variables:
                logger.info("Variable: {}".format(var))
                test_data = utils.dataset.Dataset(parameter, test=True)
                test = test_data.get_climo_variable(var, season)[:, -1]

                parameter.viewer_descr[var] = getattr(test, "long_name", var)
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = utils.general.get_name_and_yrs(
                    parameter, test_data
                )
                parameter.var_name = getattr(test, "long_name", var)
                parameter.var_units = getattr(test, "units", var)

                refs = []

                if "armdiags" in ref_name:
                    ref_file = os.path.join(
                        ref_path,
                        region[:3] + "armdiagsaciclim" + region[3:5].upper() + ".c1.nc",
                    )
                    ref_data = cdms2.open(ref_file)
                    vars_funcs = get_vars_funcs_for_derived_var(ref_data, var)
                    target_var = list(vars_funcs.keys())[0][0]
                    ref_var = ref_data(target_var)[:, 0]  # 0 mean;  1 standard devation
                    if hasattr(ref_var, "standard_name"):
                        ref_var.long_name = ref_var.standard_name
                    ref = vars_funcs[(target_var,)](utils.climo.climo(ref_var, season))

                else:
                    ref_data = utils.dataset.Dataset(parameter, ref=True)
                    ref = ref_data.get_climo_variable(var, season)[:, -1]
                ref_domain = utils.general.select_point(region, ref)
                ref.ref_name = ref_name
                refs.append(ref_domain)

                test_domain = utils.general.select_point(region, test)

                metrics_dict = create_metrics(test_domain, ref_domain)

                result = RefsTestMetrics(
                    test=test_domain, refs=refs, metrics=metrics_dict, misc=None
                )
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict["unit"] = test.units
                metrics_dict["ref_domain"] = list(ref_domain)
                metrics_dict["test_domain"] = list(test_domain)
                print(parameter.var_units, test.units)
                parameter.output_file = "-".join(
                    [ref_name, var, season, "aerosol", region]
                )
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                with open(fnm, "w") as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                fnm = os.path.join(
                    utils.general.get_output_dir(parameter.current_set, parameter),
                    parameter.output_file + ".json",
                )
                logger.info(f"Metrics saved in: {fnm}")

            if season == "ANNUALCYCLE":
                arm_diags_plot.plot_annual_cycle(var, vars_to_data[season], parameter)

    return parameter


def run_diag_pdf_daily(parameter: ARMDiagsParameter):
    logger.info("'run_diag_pdf_daily' is not yet implemented.")


def run_diag(parameter: ARMDiagsParameter) -> ARMDiagsParameter:
    if parameter.diags_set == "annual_cycle":
        return run_diag_annual_cycle(parameter)
    elif parameter.diags_set == "diurnal_cycle":
        return run_diag_diurnal_cycle(parameter)
    elif parameter.diags_set == "diurnal_cycle_zt":
        return run_diag_diurnal_cycle_zt(parameter)
    elif parameter.diags_set == "pdf_daily":
        return run_diag_pdf_daily(parameter)
    elif parameter.diags_set == "convection_onset":
        return run_diag_convection_onset(parameter)
    elif parameter.diags_set == "aerosol_activation":
        return run_diag_aerosol_activation(parameter)
    if parameter.diags_set == "annual_cycle_aerosol":
        return run_diag_annual_cycle_aerosol(parameter)
    else:
        raise Exception("Invalid diags_set={}".format(parameter.diags_set))
