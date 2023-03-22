from typing import Literal

import xarray as xr
import xcdat as xc

CLIMO_FREQ = Literal["ANN", "ANNUALCYCLE", "SEASONALCYCLE", "DJF", "MAM", "JJA", "SON"]
CDAT_TO_XCDAT_SEASON_FREQ = {
    "ANN": "month",
    "ANNUALCYCLE": "month",
    "SEASONALCYCLE": "year",
}


def climo(data_var: xr.DataArray, freq: CLIMO_FREQ) -> xr.DataArray:
    """Computes a variable's climatology for the given season.

    xCDAT's climatology API operates on a data variable within an `xr.Dataset`
    object by specifying the key of the data variable. It uses time bounds to
    redefine time as the midpoint between bounds values and month lengths for
    proper weighting.

    If the data variable is a derived variable then it is converted from an
    `xr.DataArray` to an `xr.Dataset` and time bounds are generated for the time
    axis.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.
    freq : CLIMO_FREQ
        The frequency for calculating climatology

    Returns
    -------
    xr.DataArray
        The variables' climatology
    """
    dv_key = data_var.name
    filepath = data_var.encoding.get("source")

    if filepath is not None:
        ds = xr.open_dataset(filepath)
    else:
        # The data variable is a derived variable.
        ds = data_var.to_dataset()
        ds = ds.bounds.add_bounds(axis="T")

    if freq in ["ANN", "ANNUALCYCLE", "SEASONALCYCLE"]:
        xc_freq = CDAT_TO_XCDAT_SEASON_FREQ[freq]
        ds_climo = ds.temporal.climatology(dv_key, freq=xc_freq)
    else:
        # Get the name of the time dimension and subset to the single season
        # before calculating climatology. The general best practice for
        # performance is to subset then perform calculations (split-group-apply
        # paradigm).
        time_dim = xc.get_dim_keys(data_var, axis="T")
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.season == freq)})

        ds_climo = ds.temporal.climatology(dv_key, freq="season")

    return ds_climo[dv_key]
