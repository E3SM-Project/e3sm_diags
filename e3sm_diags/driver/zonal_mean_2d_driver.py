from typing import List, Tuple

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.io import _save_data_metrics_and_plots
from e3sm_diags.driver.utils.regrid import has_z_axis, subset_and_align_datasets
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.plot.cartopy.zonal_mean_2d_plot import plot as plot_func

logger = custom_logger(__name__)


def run_diag(
    parameter: ZonalMean2dParameter, default_plevs=ZonalMean2dParameter().plevs
) -> ZonalMean2dParameter:
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = ref_ds.get_ref_climo_dataset(var_key, season, ds_test)

            # TODO: Mask ref based on reference name
            ds_ref = _mask_ref_data(ds_ref)

            # Store the variable's DataArray objects for reuse.
            dv_test = ds_test[var_key]
            dv_ref = ds_ref[var_key]

            is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
            is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

            if not is_vars_3d:
                _run_diags_2d(
                    parameter,
                    ds_test,
                    ds_ref,
                    ds_land_sea_mask,
                    season,
                    regions,
                    var_key,
                    ref_name,
                )
            elif is_vars_3d:
                _run_diags_3d(
                    parameter,
                    ds_test,
                    ds_ref,
                    ds_land_sea_mask,
                    season,
                    regions,
                    var_key,
                    ref_name,
                )
            elif is_dims_diff:
                raise RuntimeError(
                    "Dimensions of the two variables are different. Aborting."
                )

    return parameter


def _mask_ref_data(ds_ref: xr.Dataset) -> xr.Dataset:
    # TODO: Write this function
    return ds_ref


def _run_diags_2d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    for region in regions:
        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)

        (
            ds_test_region,
            ds_ref_region,
            ds_test_region_regrid,
            ds_ref_region_regrid,
            ds_diff_region,
        ) = subset_and_align_datasets(
            parameter,
            ds_test,
            ds_ref,
            ds_land_sea_mask,
            var_key,
            region,
        )

        # TODO: Mask based on variable name and reference name.
        ds_test_region_regrid, ds_ref_region_regrid = _mask_regridded_data(
            var_key, ref_name, ds_test_region_regrid, ds_ref_region_regrid
        )

        metrics_dict = _create_metrics_dict(
            var_key,
            ds_test_region,
            ds_test_region_regrid,
            ds_ref_region,
            ds_ref_region_regrid,
            ds_diff_region,
        )

        _save_data_metrics_and_plots(
            parameter,
            plot_func,
            var_key,
            ds_test_region,
            ds_ref_region,
            ds_diff_region,
            metrics_dict,
        )


def _mask_regridded_data(
    var_key: str, ref_name: str, ds_test_regrid: xr.Dataset, ds_ref_grid: xr.Dataset
) -> Tuple[xr.Dataset, xr.Dataset]:
    # TODO: Write this function
    return ds_test_regrid, ds_ref_grid


def _run_diags_3d(
    parameter: CoreParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
):
    # TODO: Write this function
    pass


def _create_metrics_dict(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_ref_regrid: xr.Dataset,
    ds_diff: xr.Dataset,
) -> MetricsDict:
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
    ds_ref : xr.Dataset
        The reference dataset.
    ds_ref_regrid : xr.Dataset
        The regridded reference dataset.
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid`` if both
        exist. This arg will be None if a model only run is performed.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    # TODO: Write this function
    metrics_dict: MetricsDict = {}

    return metrics_dict
