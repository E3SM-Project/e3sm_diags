from __future__ import annotations

import os
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd
import xarray as xr
import xcdat as xc
from scipy import interpolate

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import custom_logger
from e3sm_diags.metrics.metrics import spatial_avg
from e3sm_diags.plot import aerosol_aeronet_plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


logger = custom_logger(__name__)

# This aerosol diagnostics scripts based on AERONET sites data was originally developed by Feng Yan and adapted and integrated in e3sm_diags by Jill Zhang.
# Years include 2006–2015 average climatology for observation according to Feng et al. 2022:doi:10.1002/essoar.10510950.1, and Golaz et al. 2022 E3SMv2 paper.


def run_diag(parameter: CoreParameter) -> CoreParameter:
    """Run the aerosol aeronet diagnostics.

    Parameters
    ----------
    parameter : CoreParameter
        The parameter for the diagnostic.

    Returns
    -------
    CoreParameter
        The parameter for the diagnostic with the result (completed or failed).

    Raises
    ------
    ValueError
        If the run type is not valid.
    """
    variables = parameter.variables
    run_type = parameter.run_type
    seasons = parameter.seasons

    test_ds = Dataset(parameter, data_type="test")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            ds_test = test_ds.get_climo_dataset(var_key, season)
            da_test = ds_test[var_key]

            test_site_arr = interpolate_model_output_to_obs_sites(
                ds_test[var_key], var_key
            )

            parameter.test_name_yrs = test_ds.get_name_yrs_attr(season)
            parameter.ref_name_yrs = "AERONET (2006-2015)"

            if run_type == "model_vs_model":
                ref_ds = Dataset(parameter, data_type="ref")

                parameter.ref_name_yrs = utils.general.get_name_and_yrs(
                    parameter, ref_ds, season
                )

                ds_ref = ref_ds.get_climo_dataset(var_key, season)
                ref_site_arr = interpolate_model_output_to_obs_sites(
                    ds_ref[var_key], var_key
                )
            elif run_type == "model_vs_obs":
                ref_site_arr = interpolate_model_output_to_obs_sites(None, var_key)
            else:
                raise ValueError("Invalid run_type={}".format(run_type))

            parameter.output_file = (
                f"{parameter.ref_name}-{parameter.var_id}-{season}-global"
            )

            metrics_dict = {
                "max": da_test.max().item(),
                "min": da_test.min().item(),
                "mean": spatial_avg(ds_test, var_key, axis=["X", "Y"]),
            }
            aerosol_aeronet_plot.plot(
                parameter, da_test, test_site_arr, ref_site_arr, metrics_dict
            )

    return parameter


def interpolate_model_output_to_obs_sites(
    da_var: xr.DataArray | None, var_key: str
) -> np.ndarray:
    """Interpolate model outputs (on regular lat lon grids) to observational sites

    # TODO: Add test coverage for this function.

    Parameters
    ----------
    da_var : xr.DataArray | None
        An optional input model variable dataarray.
    var_key : str
        The key of the variable.

    Returns
    -------
    np.ndarray
        The interpolated values over all observational sites.

    Raises
    ------
    IOError
        If the variable key is invalid.
    """
    logger.info(
        "Interpolate model outputs (on regular lat lon grids) to observational sites"
    )

    if var_key == "AODABS":
        aeronet_file = os.path.join(
            e3sm_diags.INSTALL_PATH, "aerosol_aeronet/aaod550_AERONET_2006-2015.txt"
        )
        var_header = "aaod"
    elif var_key == "AODVIS":
        aeronet_file = os.path.join(
            e3sm_diags.INSTALL_PATH, "aerosol_aeronet/aod550_AERONET_2006-2015.txt"
        )
        var_header = "aod"
    else:
        raise IOError("Invalid variable input.")

    data_obs = pd.read_csv(aeronet_file, dtype=object, sep=",")

    lon_loc = np.array(data_obs["lon"].astype(float))
    lat_loc = np.array(data_obs["lat"].astype(float))
    obs_loc = np.array(data_obs[var_header].astype(float))

    num_sites = len(obs_loc)

    # Express lon_loc from 0 to 360.
    lon_loc[lon_loc < 0.0] = lon_loc[lon_loc < 0.0] + 360.0

    if da_var is not None:
        lat = xc.get_dim_coords(da_var, axis="Y")
        lon = xc.get_dim_coords(da_var, axis="X")
        f_intp = interpolate.RectBivariateSpline(lat.values, lon.values, da_var.values)

        var_intp = np.zeros(num_sites)
        for i in range(num_sites):
            var_intp[i] = f_intp(lat_loc[i], lon_loc[i])

        return var_intp

    return obs_loc
