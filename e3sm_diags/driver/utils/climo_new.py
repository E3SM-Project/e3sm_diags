from typing import Literal

import xarray as xr
import xcdat as xc

CDAT_TO_XCDAT_SEASON_FREQ = {"ANNUALCYCLE": "month", "SEASONALCYCLE": "year"}


def climo(
    ds: xr.Dataset,
    data_var: str,
    season: Literal["ANNUALCYCLE", "SEASONALCYCLE", "DJF", "MAM", "JJA", "SON"],
) -> xr.DataArray:
    """Computes a variable's climatology for the given season.

    xCDAT's climatology API uses time bounds to redefine time as the midpoint
    between bounds values and month lengths for proper weighting.

    Parameters
    ----------
    data_var : xr.Dataset
        The dataset containing the data variable and time bounds.
    data_var : str
        The name of the data variable to calculate climatology for.
    season : Literal["ANNUALCYCLE", "SEASONALCYCLE", "DJF", "MAM", "JJA", "SON"]
        The climatology season.
    """
    var_time = xc.get_dim_coords(data_var, axis="T")

    if season in ["ANNUALCYCLE", "SEASONALCYCLE"]:
        freq = CDAT_TO_XCDAT_SEASON_FREQ[season]
        ds_climo = ds.temporal.climatology(data_var, freq=freq)
    else:
        ds_climo = ds.temporal.climatology(data_var, freq="season")
        ds_climo = ds_climo.sel(f"{var_time.name}.season" == season)

    return ds_climo[data_var]
