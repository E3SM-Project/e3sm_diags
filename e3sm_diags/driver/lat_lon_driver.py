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
    align_grids_to_lower_res,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import (  # noqa: F401
    correlation,
    get_weights,
    rmse,
    spatial_avg,
    std,
)
from e3sm_diags.plot.lat_lon_plot import plot

logger = custom_logger(__name__)

# The type annotation for the metrics dictionary. The key is the
# type of metrics and the value is a sub-dictionary of metrics (key is metrics
# type and value is float). There is also a "unit" key representing the
# units for the variable.
Metrics = Dict[str, str | Dict[str, Optional[float] | List[float]]]

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Get metrics for the diagnostic.

    This function loops over each variable, season, pressure level (if 3-D),
    and region.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.

    Returns
    -------
    CoreParameter
        The parameter for the diagnostic with the result.

    Raises
    ------
    RuntimeError
        If the dimensions of the test and reference datasets are not aligned
        (e.g., one is 2-D and the other is 3-D).
    """
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    # Variables storing xarray `Dataset` objects start with `ds_` and
    # variables storing e3sm_diags `Dataset` objects end with `_ds`. This
    # is to help distinguish both objects from each other.
    test_ds = Dataset(parameter, type="test")
    ref_ds = Dataset(parameter, type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter.test_name_yrs = test_ds.get_name_and_yrs(season)
            parameter.ref_name_yrs = ref_ds.get_name_and_yrs(season)

            # The land sea mask dataset that is used for masking if the region
            # is either land or sea. This variable is instantiated here to get
            # it once per season in case it needs to be reused.
            ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

            ds_test = test_ds.get_climo_dataset(var_key, season)  # type: ignore

            # If the reference climatology dataset cannot be retrieved
            # it will be set the to the test climatology dataset which means
            # analysis is only performed on the test dataset.
            try:
                ds_ref = ref_ds.get_climo_dataset(var_key, season)  # type: ignore
                parameter.model_only = False
            except (RuntimeError, IOError):
                ds_ref = ds_test
                parameter.model_only = True

                logger.info("Cannot process reference data, analyzing test data only.")

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

            if not is_vars_3d:
                for region in regions:
                    _get_metrics_by_region(
                        parameter,
                        ds_test,
                        ds_ref,
                        var_key,
                        region,
                        ref_name,
                        season,
                        ds_land_sea_mask,
                    )
            elif is_vars_3d:
                plev = parameter.plevs
                logger.info("Selected pressure level(s): {}".format(plev))

                ds_test = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)
                ds_ref = regrid_z_axis_to_plevs(ds_ref, var_key, parameter.plevs)

                for ilev in plev:
                    ds_test_ilev = ds_test[ilev]
                    ds_ref_ilev = ds_ref[ilev]

                    for region in regions:
                        _get_metrics_by_region(
                            parameter,
                            ds_test_ilev,
                            ds_ref_ilev,
                            var_key,
                            region,
                            ref_name,
                            season,
                            ds_land_sea_mask,
                            ilev,
                        )
            elif is_dims_diff:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter


def _get_metrics_by_region(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    var_key: str,
    region: str,
    ref_name: str,
    season: str,
    ds_land_sea_mask: xr.Dataset,
    ilev: float | None = None,
):
    """Get metrics by region.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.
    ds_test : xr.Dataset
        The test dataset.
    ds_ref : xr.Dataset
        The reference dataset. If this is a model-only run then it will be the
        same dataset as ``ds_test``.
    var_key : str
        The key of the variable.
    region : str
        The region.
    ref_name : str
        The reference name.
    season : str
        The season.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is land or ocean.
    ilev : float | None, optional
        The pressure level, by default None. This option is only set if the
        variable is 3D.
    """
    logger.info(f"Selected region: {region}")
    parameter.var_region = region

    # Apply a land sea mask or subset on a specific region.
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

    # Align the grid resolutions if the diagnostic is not model only.
    if not parameter.model_only:
        ds_test_regrid, ds_ref_regrid = align_grids_to_lower_res(
            ds_test,
            ds_ref,
            var_key,
            parameter.regrid_tool,  # type: ignore
            parameter.regrid_method,
        )
        ds_diff = ds_test_regrid.copy()
        ds_diff[var_key] = ds_test_regrid[var_key] - ds_ref_regrid[var_key]
    else:
        ds_test_regrid = ds_test
        ds_ref = None
        ds_ref_regrid = None
        ds_diff = None

    metrics_dict = _create_metrics_dict(
        var_key, ds_test, ds_test_regrid, ds_ref, ds_ref_regrid, ds_diff
    )

    # Update the output and main title names based on the existence of `ilev`.
    if ilev is None:
        parameter.output_file = f"{ref_name}-{var_key}-{season}-{region}"
        parameter.main_title = f"{var_key} {season} {region}"
    else:
        ilev_str = str(int(ilev))
        parameter.output_file = f"{ref_name}-{var_key}-{ilev_str}-{season}-{region}"
        parameter.main_title = f"{var_key} {ilev_str} 'mb' {season} {region}"

    _save_and_plot_metrics(parameter, var_key, metrics_dict, ds_test, ds_ref, ds_diff)


def _create_metrics_dict(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_ref_regrid: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
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
        The regridded test Dataset. If there is no reference dataset, then this
        object is the same as ``ds_test``.
    ds_ref : xr.Dataset | None
        The optional reference dataset. This arg will be None if a model only
        run is performed.
    ds_ref_regrid : xr.Dataset | None
        The optional regridded reference dataset. This arg will be None if a
        model only run is performed.
    ds_diff : xr.Dataset | None
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid`` if both
        exist. This arg will be None if a model only run is performed.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    test_var = ds_test[var_key]
    test_regrid_var = ds_test_regrid[var_key]

    # xarray.DataArray.min() and max() returns a `np.ndarray` with a single
    # int/float element. Using `.item()` returns that single element.
    metrics_dict = {
        "test": {
            "min": test_var.min().item(),
            "max": test_regrid_var.max().item(),
            "mean": spatial_avg(ds_test, var_key, serialize=True),
        },
        "test_regrid": {
            "min": test_regrid_var.min().item(),
            "max": test_regrid_var.max().item(),
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
        ref_var = ds_ref[var_key]

        metrics_dict["ref"] = {
            "min": ref_var.min().item(),
            "max": ref_var.max().item(),
            "mean": spatial_avg(ds_ref, var_key, serialize=True),
        }

    if ds_ref_regrid is not None:
        ref_regrid_var = ds_ref_regrid[var_key]

        metrics_dict["ref_regrid"] = {
            "min": ref_regrid_var.min().item(),
            "max": ref_regrid_var.max().item(),
            "mean": spatial_avg(ds_ref_regrid, var_key, serialize=True),
            "std": std(ds_ref_regrid, var_key, serialize=True),
        }

        weights = get_weights(ds_test_regrid)
        metrics_dict["misc"] = {
            "rmse": rmse(
                test_regrid_var,
                ref_regrid_var,
                weights=weights,
                serialize=True,
            ),
            "corr": correlation(
                test_regrid_var,
                ref_regrid_var,
                weights=weights,
                serialize=True,
            ),
        }

    if ds_diff is not None:
        diff_var = ds_diff[var_key]

        metrics_dict["diff"] = {
            "min": diff_var.min().item(),
            "max": diff_var.max().item(),
            "mean": spatial_avg(ds_diff, var_key, serialize=True),
        }

    return metrics_dict


def _save_and_plot_metrics(
    parameter: CoreParameter,
    var_key: str,
    metrics_dict: Metrics,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
):
    """Save and plot the metrics.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.
    var_key : str
        The variable key.
    metrics_dict : Metrics
        The dictionary containing metrics for the variable.
    ds_test : xr.Dataset
        The test dataset.
    ds_ref : xr.Dataset | None
        The optional reference dataset. If the diagnostic is a model-only run,
        then it will be None.
    ds_diff : xr.Dataset | None
        The optional difference dataset. If the diagnostic is a model-only run,
        then it will be None.
    """
    filename = os.path.join(
        get_output_dir(parameter.current_set, parameter),
        parameter.output_file + ".json",
    )
    with open(filename, "w") as outfile:
        json.dump(metrics_dict, outfile)

    logger.info(f"Metrics saved in {filename}")

    # Set the viewer description to the "long_name" attr of the variable.
    parameter.viewer_descr[var_key] = ds_test[var_key].attrs.get(
        "long_name", "No long_name attr in test data"
    )

    plot(
        ds_test[var_key],
        ds_ref[var_key] if ds_ref is not None else None,
        ds_diff[var_key] if ds_diff is not None else None,
        metrics_dict,
        parameter,
    )

    # TODO: Write a unit test for this function call.
    if parameter.save_netcdf:
        _write_vars_to_netcdf(
            parameter,
            ds_test[var_key],
            ds_ref[var_key] if ds_ref is not None else None,
            ds_diff[var_key] if ds_diff is not None else None,
        )
