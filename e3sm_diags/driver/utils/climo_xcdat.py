"""The climatology function implemented using xCDAT's climatology API.

NOTE: This function has not been incorporated into the codebase yet because
further investigation is needed to figure out why there are large floating
point differences compared to E3SM Diags' climo function (climo.py and
climo_xr.py).
"""
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQ, CLIMO_FREQS


def climo(data_var: xr.DataArray, freq: CLIMO_FREQ) -> xr.DataArray:
    """Computes a variable's climatology for the given frequency.

    xCDAT's climatology API operates on a data variable within an `xr.Dataset`
    object by specifying the key of the data variable. It uses time bounds to
    redefine time as the midpoint between bounds values and month lengths for
    proper weighting.

    If there is no "source" filepath associated with the data variable, then it
    is considered a derived variable (a variable created from other variables
    in the dataset). We need to convert this data variable to a dataset to
    center the time coordinates and add time bounds.

    After averaging, the data variable's time dimension is squeezed and the
    time coordinates because they become singletons.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.
    freq : CLIMO_FREQ
        The frequency for calculating climatology

    Returns
    -------
    xr.DataArray
        The variable's climatology
    """
    dv_key = data_var.name
    filepath = data_var.encoding.get("source")
    time_dim = xc.get_dim_keys(data_var, axis="T")

    if filepath is None:
        ds = data_var.to_dataset()
        ds = xc.center_times(ds)
        ds = ds.bounds.add_bounds(axis="T")
    else:
        ds = xc.open_dataset(filepath, center_times=True)

    if freq in ["ANN"]:
        ds_climo = ds.temporal.average(dv_key, weighted=True)
    elif freq in [
        "01",
        "02",
        "03",
        "04",
        "05",
        "06",
        "07",
        "08",
        "09",
        "10",
        "11",
        "12",
    ]:
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.month == freq)})
        ds_climo = ds.temporal.climatology(dv_key, freq="month", weighted=True)
    elif freq in ["DJF", "MAM", "JJA", "SON"]:
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.season == freq)})
        ds_climo = ds.temporal.climatology(dv_key, freq="season", weighted=True)
    else:
        raise ValueError(
            f"`freq='{freq}'` is not a valid climatology frequency. Options "
            f"include {CLIMO_FREQS}'"
        )

    dv_climo = ds_climo[dv_key].copy()
    dv_climo = dv_climo.squeeze(dim=time_dim)
    dv_climo = dv_climo.drop_vars(time_dim)

    return dv_climo
