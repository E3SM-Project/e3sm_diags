from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING, Dict, Optional

import xarray as xr

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.driver.utils.io import _write_vars_to_netcdf
from e3sm_diags.driver.utils.regrid import convert_z_axis_to_pressure_levels, has_z_axis
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std  # noqa: F401
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
        for var_key in variables:
            logger.info("Variable: {}".format(var_key))
            parameter.var_id = var_key

            ds_climo_test = test_ds.get_climo_dataset(var_key, season)  # type: ignore

            try:
                ds_climo_ref = ref_ds.get_climo_dataset(var_key, season)  # type: ignore
            except (RuntimeError, IOError):
                ds_climo_ref = ds_climo_test

                logger.info("Cannot process reference data, analyzing test data only.")
                parameter.model_only = True

            # Store the variable's DataArray objects for reuse.
            dv_climo_test = ds_climo_test[var_key]
            dv_climo_ref = ds_climo_ref[var_key]

            # Set the viewer description.
            parameter.viewer_descr[var_key] = dv_climo_ref.attrs.get(
                "long_name", "No long_name attr in test data"
            )

            # Check whether the variable has a Z axis. If it does not, just
            # save metrics for each region to a JSON file. If it does,
            # then convert the Z axis to pressure levels before saving metrics
            # for each region to a JSON file.
            vars_have_z_axis = has_z_axis(dv_climo_test) and has_z_axis(dv_climo_ref)

            if not vars_have_z_axis:
                for region in regions:
                    logger.info(f"Selected region: {region}")

                    parameter.var_region = region
                    parameter.output_file = f"{ref_name}-{var_key}-{season}-{region}"
                    parameter.main_title = f"{var_key} {season} {region}"

                    mv1_domain = utils.general.select_region(
                        region, dv_climo_test, land_frac, ocean_frac, parameter
                    )
                    mv2_domain = utils.general.select_region(
                        region, dv_climo_ref, land_frac, ocean_frac, parameter
                    )

                    create_and_save_data_and_metrics(parameter, mv1_domain, mv2_domain)
            elif vars_have_z_axis:
                plev = parameter.plevs
                logger.info("Selected pressure level: {}".format(plev))

                mv1_p = convert_z_axis_to_pressure_levels(
                    ds_climo_test, var_key, parameter.plevs
                )
                mv2_p = convert_z_axis_to_pressure_levels(
                    ds_climo_ref, var_key, parameter.plevs
                )

                # Loop over each pressure level, subset the variable on that
                # pressure level, subset on the region, then save the related
                # metrics.
                for ilev in range(len(plev)):
                    dv_climo_test = mv1_p[ilev]
                    dv_climo_ref = mv2_p[ilev]

                    for region in regions:
                        logger.info(f"Selected region: {region}")

                        parameter.var_region = region
                        parameter.output_file = f"{ref_name}-{var_key}-{str(int(plev[ilev]))}-{season}-{region}"
                        parameter.main_title = (
                            f"{var_key} {str(int(plev[ilev]))} 'mb' {season} {region}"
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
            elif has_z_axis(dv_climo_test) != has_z_axis(dv_climo_ref):
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

        reg_diff = mv1_reg - mv2_reg
    else:
        mv1_reg = mv1_domain
        mv2_domain = None
        mv2_reg = None
        reg_diff = None

    metrics_dict = _create_metrics(mv1_domain, mv1_reg, mv2_domain, mv2_reg, reg_diff)

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
        reg_diff,
        metrics_dict,
        parameter,
    )

    # TODO: Write a unit test for this function call.
    _write_vars_to_netcdf(
        parameter,
        mv1_domain,
        mv2_domain,
        reg_diff,
    )


def _create_metrics(
    test: xr.DataArray,
    test_regrid: xr.DataArray,
    ref: Optional[xr.DataArray],
    ref_regrid: Optional[xr.DataArray],
    diff_regrid: Optional[xr.DataArray],
) -> Dict[str, Dict]:
    """Computes metrics for variables.

    Metrics include min value, max value, spatial average (mean), standard
    deviation, correlation (pearson_r), and RMSE. The default value for
    optional arguments is 999.999, which represents missing metrics.

    Parameters
    ----------
    test : xr.DataArray
        The test variable, which is usually a "domain" variable.
    test_regrid : xr.DataArray
        The regridded test variable.
    ref : Optional[xr.DataArray]
        The optional reference variable, which is usually a "domain" variable.
    ref_regrid : Optional[xr.DataArray]
        The optional regridded reference variable.
    diff_regrid : Optional[xr.DataArray]
        The optional difference between the regridded variables.

    Returns
    -------
    Dict[str, Dict]
        A dictionary with the key being the name of the parameter (and "misc")
        and the value being the related metrics. Example:
        {
            "test": {
                "min": test.min().item(),
                "max": test.max().item(),
            }
        }
    """
    # TODO: Figure out how to handle arguments for this function. Do we need
    # to pass dataset objects for spatial avg and standard deviation?
    default_value = 999.999

    # xarray.DataArray.min() and max() returns a `np.ndarray` with a single
    # int/float element. Using `.item()` returns that single element.
    metrics_dict: Dict[str, Dict] = {
        "test": {
            "min": test.min().item(),
            "max": test.max().item(),
            # "mean": spatial_avg(ds_test, test, serialize=True),
        },
        "test_regrid": {
            "min": test_regrid.min().item(),
            "max": test_regrid.max().item(),
            # "mean": spatial_avg(ds_test, test_regrid, serialize=True),
            # "std": std(ds_test, test_regrid, serialize=True),
        },
        "ref": {
            "min": default_value,
            "max": default_value,
            "mean": default_value,
        },
        "ref_regrid": {
            "min": default_value,
            "max": default_value,
            "mean": default_value,
            "std": default_value,
        },
        "misc": {
            "rmse": default_value,
            "corr": default_value,
        },
    }

    if ref is not None:
        metrics_dict["ref"] = {
            "min": ref.min().item(),
            "max": ref.max().item(),
            # "mean": spatial_avg(ds_ref, ref, serialize=True),
        }

    if ref_regrid is not None:
        metrics_dict["ref_regrid"] = {
            "min": ref_regrid.min().item(),
            "max": ref_regrid.max().item(),
            # "mean": spatial_avg(ref_regrid, serialize=True),
            # "std": std(ref_regrid, serialize=True),
        }
        metrics_dict["misc"] = {
            "rmse": rmse(test_regrid, ref_regrid, serialize=True),
            "corr": correlation(test_regrid, ref_regrid, serialize=True),
        }

    if diff_regrid is not None:
        metrics_dict["diff"] = {
            "min": diff_regrid.min().item(),
            "max": diff_regrid.max().item(),
            # "mean": spatial_avg(diff_regrid, serialize=True),
        }

    return metrics_dict
