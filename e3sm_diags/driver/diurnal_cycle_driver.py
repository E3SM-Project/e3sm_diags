from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import xarray as xr

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.diurnal_cycle_xr import composite_diurnal_cycle
from e3sm_diags.driver.utils.io import _get_output_filename_filepath
from e3sm_diags.driver.utils.regrid import _apply_land_sea_mask, _subset_on_region
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.diurnal_cycle_plot import plot as plot_func

logger = custom_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter


def run_diag(parameter: DiurnalCycleParameter) -> DiurnalCycleParameter:
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
            ds_ref = ref_ds.get_climo_dataset(var_key, season)

            for region in regions:
                if "land" in region or "ocean" in region:
                    test_domain = _apply_land_sea_mask(
                        ds_test,
                        ds_land_sea_mask,
                        var_key,
                        region,  # type: ignore
                        parameter.regrid_tool,
                        parameter.regrid_method,
                    )

                    ref_domain = _apply_land_sea_mask(
                        ds_ref,
                        ds_land_sea_mask,
                        var_key,
                        region,  # type: ignore
                        parameter.regrid_tool,
                        parameter.regrid_method,
                    )
                else:
                    test_domain = ds_test.copy()
                    ref_domain = ds_ref.copy()

                test_domain = _subset_on_region(test_domain, var_key, region)
                ref_domain = _subset_on_region(ref_domain, var_key, region)

                parameter.viewer_descr[var_key] = ds_test[var_key].attrs.get(
                    "long_name", "No long_name attr in test data."
                )
                parameter.output_file = "-".join([ref_name, var_key, season, region])
                parameter.main_title = str(
                    " ".join([var_key, "Diurnal Cycle ", season, region])
                )

                (
                    test_cmean,
                    test_amplitude,
                    test_maxtime,
                ) = composite_diurnal_cycle(
                    test_domain, var_key, season
                )  # type: ignore
                (
                    ref_cmean,
                    ref_amplitude,
                    ref_maxtime,
                ) = composite_diurnal_cycle(
                    ref_domain, var_key, season
                )  # type: ignore

                parameter.var_region = region

                plot_func(
                    test_maxtime,
                    test_amplitude,
                    ref_maxtime,
                    ref_amplitude,
                    parameter,
                )

                ds_out_test = xr.Dataset(
                    data_vars={
                        test_cmean.name: test_cmean,
                        test_amplitude.name: test_amplitude,
                        test_maxtime.name: test_maxtime,
                    }
                )
                ds_out_ref = xr.Dataset(
                    data_vars={
                        ref_cmean.name: ref_cmean,
                        ref_amplitude.name: ref_amplitude,
                        ref_maxtime.name: ref_maxtime,
                    }
                )

                _write_vars_to_netcdf(parameter, var_key, ds_out_test, "test")
                _write_vars_to_netcdf(parameter, var_key, ds_out_ref, "ref")

    return parameter


def _write_vars_to_netcdf(
    parameter: DiurnalCycleParameter,
    var_key: str,
    ds: xr.Dataset,
    data_type: Literal["test", "ref"],
):
    _, filepath = _get_output_filename_filepath(parameter, data_type)

    ds.to_netcdf(filepath)

    logger.info(f"'{var_key}' {data_type} variable output saved in: {filepath}")
