from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING

import xarray as xr

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.driver.utils.regrid import (
    convert_to_pressure_levels,
    has_z_axis_coords,
)
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics import corr, mean, rmse, std
from e3sm_diags.plot import plot

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:  # noqa: C901
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    test_ds = Dataset(parameter, type="test")
    ref_ds = Dataset(parameter, type="ref")

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = test_ds.get_name_and_yrs(season)
        parameter.ref_name_yrs = ref_ds.get_name_and_yrs(season)

        try:
            ds_land_frac = test_ds.get_climo_dataset("LANDFRAC", season)  # type: ignore
            land_frac = ds_land_frac["LANDFRAC"]

            ds_ocean_frac = test_ds.get_climo_dataset("OCNFRAC", season)  # type: ignore
            ocean_frac = ds_ocean_frac["OCNFRAC"]
        except RuntimeError as e:
            logger.warning(e)

            mask_path = os.path.join(
                e3sm_diags.INSTALL_PATH, "acme_ne30_ocean_land_mask.nc"
            )

            ds_mask = xr.open_dataset(mask_path)
            land_frac = ds_mask["LANDFRAC"]
            ocean_frac = ds_mask["OCNFRAC"]

        parameter.model_only = False
        for var in variables:
            logger.info("Variable: {}".format(var))
            parameter.var_id = var

            ds_climo_test = test_ds.get_climo_dataset(var, season)  # type: ignore
            dv_climo_test = ds_climo_test[var]

            try:
                ds_climo_ref = ref_ds.get_climo_dataset(var, season)  # type: ignore
                dv_climo_ref = ds_climo_ref[var]
            except (RuntimeError, IOError):
                ds_climo_ref = ds_climo_test
                dv_climo_ref = dv_climo_test.copy()
                logger.info("Can not process reference data, analyse test data only")

                parameter.model_only = True

            # Set the viewer description.
            parameter.viewer_descr[var] = dv_climo_test.attrs.get(
                "long_name", "No long_name attr in test data"
            )

            # For variables with a z-axis.
            if has_z_axis_coords(dv_climo_test) and has_z_axis_coords(dv_climo_ref):
                plev = parameter.plevs
                logger.info("Selected pressure level: {}".format(plev))

                mv1_p = convert_to_pressure_levels(
                    ds_climo_test, dv_climo_test, parameter.plevs
                )
                mv2_p = convert_to_pressure_levels(
                    ds_climo_ref, dv_climo_ref, parameter.plevs
                )

                # Select plev.
                for ilev in range(len(plev)):
                    dv_climo_test = mv1_p[ilev]
                    dv_climo_ref = mv2_p[ilev]

                    for region in regions:
                        logger.info(f"Selected region: {region}")

                        parameter.var_region = region
                        parameter.output_file = (
                            f"{ref_name}-{var}-{str(int(plev[ilev]))}-{season}-{region}"
                        )
                        parameter.main_title = (
                            f"{var} {str(int(plev[ilev]))} 'mb' {season} {region}"
                        )

                        mv1_domain = utils.general.select_region(
                            region, dv_climo_test, land_frac, ocean_frac, parameter
                        )
                        mv2_domain = utils.general.select_region(
                            region, dv_climo_ref, land_frac, ocean_frac, parameter
                        )

                        create_and_save_data_and_metrics(
                            parameter, mv1_domain, mv2_domain
                        )

            # For variables without a z-axis.
            elif not has_z_axis_coords(dv_climo_test) and not has_z_axis_coords(
                dv_climo_ref
            ):
                for region in regions:
                    logger.info(f"Selected region: {region}")

                    parameter.var_region = region
                    parameter.output_file = f"{ref_name}-{var}-{season}-{region}"
                    parameter.main_title = f"{var} {season} {region}"

                    mv1_domain = utils.general.select_region(
                        region, dv_climo_test, land_frac, ocean_frac, parameter
                    )
                    mv2_domain = utils.general.select_region(
                        region, dv_climo_ref, land_frac, ocean_frac, parameter
                    )

                    create_and_save_data_and_metrics(parameter, mv1_domain, mv2_domain)
            else:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter


def create_and_save_data_and_metrics(parameter, mv1_domain, mv2_domain):
    if not parameter.model_only:
        # Regrid towards the lower resolution of the two
        # variables for calculating the difference.
        mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
            mv1_domain,
            mv2_domain,
            parameter.regrid_tool,
            parameter.regrid_method,
        )

        diff = mv1_reg - mv2_reg
    else:
        mv2_domain = None
        mv2_reg = None
        mv1_reg = mv1_domain
        diff = None

    metrics_dict = create_metrics(mv2_domain, mv1_domain, mv2_reg, mv1_reg, diff)

    # Saving the metrics as a json.
    metrics_dict["unit"] = mv1_domain.units

    fnm = os.path.join(
        utils.general.get_output_dir(parameter.current_set, parameter),
        parameter.output_file + ".json",
    )
    with open(fnm, "w") as outfile:
        json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in {fnm}")

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


def create_metrics(
    ref: xr.DataArray,
    test: xr.DataArray,
    ref_regrid: xr.DataArray,
    test_regrid: xr.DataArray,
    diff: xr.DataArray,
):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    # For input None, metrics are instantiated to 999.999.
    # Apply float() to make sure the elements in metrics_dict are JSON serializable, i.e. np.float64 type is JSON serializable, but not np.float32.
    missing_value = 999.999
    metrics_dict = {}
    metrics_dict["ref"] = {
        "min": float(ref.min()) if ref is not None else missing_value,
        "max": float(ref.max()) if ref is not None else missing_value,
        "mean": float(mean(ref)) if ref is not None else missing_value,
    }
    metrics_dict["ref_regrid"] = {
        "min": float(ref_regrid.min()) if ref_regrid is not None else missing_value,
        "max": float(ref_regrid.max()) if ref_regrid is not None else missing_value,
        "mean": float(mean(ref_regrid)) if ref_regrid is not None else missing_value,
        "std": float(std(ref_regrid)) if ref_regrid is not None else missing_value,
    }
    metrics_dict["test"] = {
        "min": float(test.min()),
        "max": float(test.max()),
        "mean": float(mean(test)),
    }
    metrics_dict["test_regrid"] = {
        "min": float(test_regrid.min()),
        "max": float(test_regrid.max()),
        "mean": float(mean(test_regrid)),
        "std": float(std(test_regrid)),
    }
    metrics_dict["diff"] = {
        "min": float(diff.min()) if diff is not None else missing_value,
        "max": float(diff.max()) if diff is not None else missing_value,
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
