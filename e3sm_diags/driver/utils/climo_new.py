from typing import Literal, get_args

import xarray as xr
import xcdat as xc

CLIMO_FREQ = Literal["ANN", "ANNUALCYCLE", "SEASONALCYCLE", "DJF", "MAM", "JJA", "SON"]


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

    # If there is no filepath associated with the data variable, then it is
    # considered a derived variable. We need to convert this data variable to a
    # dataset to center the time coordinates and add time bounds.
    if filepath is None:
        ds = data_var.to_dataset()
        ds = xc.center_times(ds)
        ds = ds.bounds.add_bounds(axis="T")
    else:
        ds = xc.open_dataset(filepath, center_times=True)

    if freq in ["ANNUALCYCLE", "ANN"]:
        ds_climo = ds.temporal.average(dv_key, weighted=True)
    elif freq == "SEASONALCYCLE":
        ds_climo = ds.temporal.climatology(dv_key, freq="season", weighted=True)
    elif freq in ["DJF", "MAM", "JJA", "SON"]:
        # The general best practice for performance is to subset then perform
        # calculations (split-group-apply paradigm).
        time_dim = xc.get_dim_keys(data_var, axis="T")
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.season == freq)})

        ds_climo = ds.temporal.climatology(dv_key, freq="season", weighted=True)
    else:
        raise ValueError(
            f"`freq='{freq}'` is not a valid climatology frequency. Options "
            f"include {get_args(CLIMO_FREQ)}'"
        )

    return ds_climo[dv_key]
