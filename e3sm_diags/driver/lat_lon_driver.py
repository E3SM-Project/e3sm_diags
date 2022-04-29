from __future__ import print_function

import json
import os

import cdms2

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics import corr, max_cdms, mean, min_cdms, rmse, std
from e3sm_diags.plot import plot

logger = custom_logger(__name__)


def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    missing_value = 999.999
    metrics_dict = {}
    metrics_dict["ref"] = {
        "min": float(min_cdms(ref)) if ref is not None else missing_value,
        "max": float(max_cdms(ref)) if ref is not None else missing_value,
        "mean": float(mean(ref)) if ref is not None else missing_value,
    }
    metrics_dict["ref_regrid"] = {
        "min": float(min_cdms(ref_regrid)) if ref_regrid is not None else missing_value,
        "max": float(max_cdms(ref_regrid)) if ref_regrid is not None else missing_value,
        "mean": float(mean(ref_regrid)) if ref_regrid is not None else missing_value,
        "std": float(std(ref_regrid)) if ref_regrid is not None else missing_value,
    }
    metrics_dict["test"] = {
        "min": float(min_cdms(test)),
        "max": float(max_cdms(test)),
        "mean": float(mean(test)),
    }
    metrics_dict["test_regrid"] = {
        "min": float(min_cdms(test_regrid)),
        "max": float(max_cdms(test_regrid)),
        "mean": float(mean(test_regrid)),
        "std": float(std(test_regrid)),
    }
    metrics_dict["diff"] = {
        "min": float(min_cdms(diff)) if diff is not None else missing_value,
        "max": float(max_cdms(diff)) if diff is not None else missing_value,
        "mean": float(mean(diff)) if diff is not None else missing_value,
    }
    metrics_dict["misc"] = {
        "rmse": float(rmse(test_regrid, ref_regrid))
        if ref_regrid is not None
        else missing_value,
        "corr": float(corr(test_regrid, ref_regrid))
        if ref_regrid is not None
        else missing_value,
    }
    return metrics_dict


def run_diag(parameter):  # noqa: C901
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(
            parameter, test_data, season
        )
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(
            parameter, ref_data, season
        )

        # Get land/ocean fraction for masking.
        try:
            land_frac = test_data.get_climo_variable("LANDFRAC", season)
            ocean_frac = test_data.get_climo_variable("OCNFRAC", season)
        except Exception:
            mask_path = os.path.join(
                e3sm_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc"
            )
            with cdms2.open(mask_path) as f:
                land_frac = f("LANDFRAC")
                ocean_frac = f("OCNFRAC")

        for var in variables:
            logger.info("Variable: {}".format(var))
            parameter.var_id = var

            mv1 = test_data.get_climo_variable(var, season)
            try:
                mv2 = ref_data.get_climo_variable(var, season)
            except Exception:
                mv2 = mv1
                logger.info("Can not process reference data, analyse test data only")
                parameter.ref_name = ""
                parameter.case_id = "Model_only"

            parameter.viewer_descr[var] = (
                mv1.long_name
                if hasattr(mv1, "long_name")
                else "No long_name attr in test data."
            )

            # For variables with a z-axis.
            if mv1.getLevel() and mv2.getLevel():
                plev = parameter.plevs
                logger.info("Selected pressure level: {}".format(plev))

                mv1_p = utils.general.convert_to_pressure_levels(
                    mv1, plev, test_data, var, season
                )
                mv2_p = utils.general.convert_to_pressure_levels(
                    mv2, plev, ref_data, var, season
                )

                # Select plev.
                for ilev in range(len(plev)):
                    mv1 = mv1_p[
                        ilev,
                    ]
                    mv2 = mv2_p[
                        ilev,
                    ]

                    for region in regions:
                        logger.info(f"Selected regions: {region}")
                        mv1_domain = utils.general.select_region(
                            region, mv1, land_frac, ocean_frac, parameter
                        )
                        mv2_domain = utils.general.select_region(
                            region, mv2, land_frac, ocean_frac, parameter
                        )

                        parameter.output_file = "-".join(
                            [
                                ref_name,
                                var,
                                str(int(plev[ilev])),
                                season,
                                region,
                            ]
                        )
                        parameter.main_title = str(
                            " ".join(
                                [
                                    var,
                                    str(int(plev[ilev])),
                                    "mb",
                                    season,
                                    region,
                                ]
                            )
                        )

                        if parameter.ref_name != "":
                            # Regrid towards the lower resolution of the two
                            # variables for calculating the difference.
                            mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                                mv1_domain,
                                mv2_domain,
                                parameter.regrid_tool,
                                parameter.regrid_method,
                            )

                            diff = mv1_reg - mv2_reg
                            metrics_dict = create_metrics(
                                mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff
                            )
                        else:
                            mv2_domain = None
                            mv2_reg = None
                            diff = None
                            metrics_dict = create_metrics(
                                mv2_domain, mv1, mv2_reg, mv1, diff
                            )

                            # Saving the metrics as a json.
                        metrics_dict["unit"] = mv1.units
                        fnm = os.path.join(
                            utils.general.get_output_dir(
                                parameter.current_set, parameter
                            ),
                            parameter.output_file + ".json",
                        )
                        with open(fnm, "w") as outfile:
                            json.dump(metrics_dict, outfile)
                        # Get the filename that the user has passed in and display that.
                        fnm = os.path.join(
                            utils.general.get_output_dir(
                                parameter.current_set, parameter
                            ),
                            parameter.output_file + ".json",
                        )
                        logger.info(f"Metrics saved in: {fnm}")

                        parameter.var_region = region
                        plot(
                            parameter.current_set,
                            mv2_domain,
                            mv1_domain,
                            diff,
                            metrics_dict,
                            parameter,
                        )
                        utils.general.save_ncfiles(
                            parameter.current_set,
                            mv1_domain,
                            mv2_domain,
                            diff,
                            parameter,
                        )

            # For variables without a z-axis.
            elif mv1.getLevel() is None and mv2.getLevel() is None:
                for region in regions:
                    logger.info(f"Selected region: {region}")
                    mv1_domain = utils.general.select_region(
                        region, mv1, land_frac, ocean_frac, parameter
                    )
                    mv2_domain = utils.general.select_region(
                        region, mv2, land_frac, ocean_frac, parameter
                    )

                    parameter.output_file = "-".join([ref_name, var, season, region])
                    parameter.main_title = str(" ".join([var, season, region]))

                    if parameter.ref_name != "":
                        # Regrid towards the lower resolution of the two
                        # variables for calculating the difference.
                        mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                            mv1_domain,
                            mv2_domain,
                            parameter.regrid_tool,
                            parameter.regrid_method,
                        )

                        diff = mv1_reg - mv2_reg
                        metrics_dict = create_metrics(
                            mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff
                        )
                    else:
                        mv2_domain = None
                        mv2_reg = None
                        diff = None
                        metrics_dict = create_metrics(
                            mv2_domain, mv1, mv2_reg, mv1, diff
                        )

                    # Saving the metrics as a json.
                    metrics_dict["unit"] = mv1.units
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
                    logger.info(f"Metrics saved in {fnm}")

                    parameter.var_region = region
                    plot(
                        parameter.current_set,
                        mv2_domain,
                        mv1_domain,
                        diff,
                        metrics_dict,
                        parameter,
                    )
                    utils.general.save_ncfiles(
                        parameter.current_set,
                        mv1_domain,
                        mv2_domain,
                        diff,
                        parameter,
                    )

            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter
