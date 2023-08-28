from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING, Dict, List, Optional

import xarray as xr

from e3sm_diags.driver.utils.dataset_new import Dataset
from e3sm_diags.driver.utils.general import get_output_dir
from e3sm_diags.driver.utils.io import _write_vars_to_netcdf
from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
    has_z_axis,
    regrid_to_lower_res,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import correlation, rmse, spatial_avg, std  # noqa: F401
from e3sm_diags.plot.lat_lon_plot import plot

logger = custom_logger(__name__)

# The type annotation for the metrics dictionary. The key is the
# type of metrics and the value is a sub-dictionary of metrics (key is metrics
# type and value is float). There is also a "unit" key representing the
# units for the variable.
Metrics = Dict[str, str | Dict[str, Optional[float] | List[float]]]

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:  # noqa: C901
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    # NOTE: There is a naming conflict between the e3sm_diags `Dataset` class
    # and the xarray `Dataset` class, which makes referencing both classes
    # potentially confusing.
    test_ds = Dataset(parameter, type="test")
    ref_ds = Dataset(parameter, type="ref")

    # TODO: Add notes on how this for loop works.
    # Loop over each season and variable to calculate metrics. If a variable
    # has a Z axis, it is regridded to pressure levels then loop over
    # the pressure levels to calculate the metrics.

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = test_ds.get_name_and_yrs(season)
        parameter.ref_name_yrs = ref_ds.get_name_and_yrs(season)

        # The land sea mask used during regional selection.
        ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

        parameter.model_only = False
        for var_key in variables:
            logger.info("Variable: {}".format(var_key))
            parameter.var_id = var_key

            ds_test = test_ds.get_climo_dataset(var_key, season)  # type: ignore

            try:
                ds_ref = ref_ds.get_climo_dataset(var_key, season)  # type: ignore
            except (RuntimeError, IOError):
                ds_ref = ds_test

                logger.info("Cannot process reference data, analyzing test data only.")
                parameter.model_only = True

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            # Set the viewer description.
            parameter.viewer_descr[var_key] = dv_test.attrs.get(
                "long_name", "No long_name attr in test data"
            )

            # Check whether the variable has a Z axis. If it does not, just
            # save metrics for each region to a JSON file. If it does,
            # then convert the Z axis to pressure levels before saving metrics
            # for each region to a JSON file.
            vars_have_z_axis = has_z_axis(dv_test) and has_z_axis(dv_ref)

            # TODO: Refactor both conditionals since logic is similar.
            if not vars_have_z_axis:
                for region in regions:
                    logger.info(f"Selected region: {region}")

                    parameter.var_region = region
                    parameter.output_file = f"{ref_name}-{var_key}-{season}-{region}"
                    parameter.main_title = f"{var_key} {season} {region}"

                    if region == "land" or region == "ocean":
                        ds_test = _apply_land_sea_mask(
                            ds_test,
                            ds_land_sea_mask,
                            var_key,
                            region,  # type: ignore
                            parameter.regrid_tool,
                            parameter.regrid_method,
                        )
                        ds_ref = _apply_land_sea_mask(
                            ds_ref,
                            ds_land_sea_mask,
                            var_key,
                            region,  # type: ignore
                            parameter.regrid_tool,
                            parameter.regrid_method,
                        )
                    elif region != "global":
                        ds_test = _subset_on_region(ds_test, var_key, region)
                        ds_ref = _subset_on_region(ds_ref, var_key, region)

                    create_and_save_data_and_metrics(
                        parameter, ds_test, ds_ref, var_key
                    )
            elif vars_have_z_axis:
                plev = parameter.plevs
                logger.info("Selected pressure level: {}".format(plev))

                mv1_p = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)
                mv2_p = regrid_z_axis_to_plevs(ds_ref, var_key, parameter.plevs)

                # Loop over each pressure level, subset the variable on that
                # pressure level, subset on the region, then save the related
                # metrics.
                for ilev in range(len(plev)):
                    dv_test = mv1_p[ilev]
                    dv_ref = mv2_p[ilev]

                    for region in regions:
                        logger.info(f"Selected region: {region}")

                        parameter.var_region = region
                        parameter.output_file = f"{ref_name}-{var_key}-{str(int(plev[ilev]))}-{season}-{region}"
                        parameter.main_title = (
                            f"{var_key} {str(int(plev[ilev]))} 'mb' {season} {region}"
                        )

                        if region == "land" or region == "ocean":
                            ds_test = _apply_land_sea_mask(
                                ds_test,
                                ds_land_sea_mask,
                                var_key,
                                region,  # type: ignore
                                parameter.regrid_tool,
                                parameter.regrid_method,
                            )
                            ds_ref = _apply_land_sea_mask(
                                ds_ref,
                                ds_land_sea_mask,
                                var_key,
                                region,  # type: ignore
                                parameter.regrid_tool,
                                parameter.regrid_method,
                            )
                        elif region != "global":
                            ds_test = _subset_on_region(ds_test, var_key, region)
                            ds_ref = _subset_on_region(ds_ref, var_key, region)

                        create_and_save_data_and_metrics(
                            parameter, ds_test, ds_ref, var_key
                        )
            elif has_z_axis(dv_test) != has_z_axis(dv_ref):
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter


def create_and_save_data_and_metrics(
    parameter: CoreParameter, ds_test: xr.Dataset, ds_ref: xr.Dataset, var_key: str
):
    if parameter.model_only:
        ds_test_regrid = ds_ref
        ds_ref = None
        ds_ref_regrid = None
        ds_diff = None
    else:
        # To calculate the difference between the datasets, their grids are
        # aligned the grid of the dataset with the lower resolution. No regridding
        # is performed if both datasets have the same resolution.
        ds_test_regrid, ds_ref_regrid = regrid_to_lower_res(
            ds_test,
            ds_ref,
            var_key,
            parameter.regrid_tool,  # type: ignore
            parameter.regrid_method,
        )

        ds_diff = ds_test_regrid.copy()
        ds_diff[var_key] = ds_test_regrid[var_key] - ds_ref_regrid[var_key]

    metrics_dict = _create_metrics(
        var_key, ds_test, ds_test_regrid, ds_ref, ds_ref_regrid, ds_diff
    )

    # Saving the metrics as a json.
    filename = os.path.join(
        get_output_dir(parameter.current_set, parameter),
        parameter.output_file + ".json",
    )
    with open(filename, "w") as outfile:
        json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in {filename}")

    plot(
        ds_test[var_key],
        ds_ref[var_key] if ds_ref is not None else None,
        ds_diff[var_key] if ds_diff is not None else None,
        metrics_dict,
        parameter,
    )

    # TODO: Write a unit test for this function call.
    _write_vars_to_netcdf(
        parameter,
        ds_test[var_key],
        ds_ref[var_key] if ds_ref is not None else None,
        ds_diff[var_key] if ds_diff is not None else None,
    )


def _create_metrics(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset,
    ds_ref: Optional[xr.Dataset],
    ds_ref_regrid: Optional[xr.Dataset],
    ds_diff: Optional[xr.Dataset],
) -> Metrics:
    """Calculate metrics using the variable in the datasets.

    Metrics include min value, max value, spatial average (mean), standard
    deviation, correlation (pearson_r), and RMSE. The default value for
    optional metrics is None.

    Parameters
    ----------
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_test_regrid : xr.Dataset
        The regridded test Dataset.
    ds_ref : Optional[xr.Dataset]
        The optional reference dataset.
    ds_ref_regrid : Optional[xr.Dataset]
        The optional regridded reference dataset.
    ds_diff : Optional[xr.Dataset]
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid``.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    # xarray.DataArray.min() and max() returns a `np.ndarray` with a single
    # int/float element. Using `.item()` returns that single element.
    metrics_dict: Metrics = {
        "test": {
            "min": ds_test[var_key].min().item(),
            "max": ds_test[var_key].max().item(),
            "mean": spatial_avg(ds_test, var_key, serialize=True),
        },
        "test_regrid": {
            "min": ds_test_regrid[var_key].min().item(),
            "max": ds_test_regrid[var_key].max().item(),
            "mean": spatial_avg(ds_test_regrid, var_key, serialize=True),
            "std": std(ds_test_regrid, var_key, serialize=True),
        },
        "ref": {
            "min": None,
            "max": None,
            "mean": None,
        },
        "ref_regrid": {
            "min": None,
            "max": None,
            "mean": None,
            "std": None,
        },
        "misc": {
            "rmse": None,
            "corr": None,
        },
        "diff": {
            "min": None,
            "max": None,
            "mean": None,
        },
        "unit": ds_test[var_key].attrs["units"],
    }

    if ds_ref is not None:
        metrics_dict["ref"] = {
            "min": ds_ref[var_key].min().item(),
            "max": ds_ref[var_key].max().item(),
            "mean": spatial_avg(ds_ref, var_key, serialize=True),
        }

    if ds_ref_regrid is not None:
        metrics_dict["ref_regrid"] = {
            "min": ds_ref_regrid[var_key].min().item(),
            "max": ds_ref_regrid[var_key].max().item(),
            "mean": spatial_avg(ds_ref_regrid, var_key, serialize=True),
            "std": std(ds_ref_regrid, var_key, serialize=True),
        }
        metrics_dict["misc"] = {
            "rmse": rmse(
                ds_test_regrid[var_key], ds_ref_regrid[var_key], serialize=True
            ),
            "corr": correlation(
                ds_test_regrid[var_key], ds_ref_regrid[var_key], serialize=True
            ),
        }

    if ds_diff is not None:
        metrics_dict["diff"] = {
            "min": ds_diff[var_key].min().item(),
            "max": ds_diff[var_key].max().item(),
            "mean": spatial_avg(ds_diff, var_key, serialize=True),
        }

    return metrics_dict
