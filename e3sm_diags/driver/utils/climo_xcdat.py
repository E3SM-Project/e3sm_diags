"""The climatology function implemented using xCDAT's climatology API.

NOTE: This function has not been incorporated into the codebase yet because
further investigation is needed to figure out why there are large floating
point differences compared to E3SM Diags' climo function (climo.py and
climo_xr.py).
"""
import xarray as xr
import xcdat as xc

from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQ, CLIMO_FREQS


def climo(dataset: xr.Dataset, var_key: str, freq: CLIMO_FREQ) -> xr.DataArray:
    """Computes a variable's climatology for the given frequency.

    xCDAT's climatology API operates on a data variable within an `xr.Dataset`
    object by specifying the key of the data variable. It uses time bounds to
    redefine time as the midpoint between bounds values and month lengths for
    proper weighting.

    After averaging, the data variable's time dimension is squeezed and the
    time coordinates because they become singletons.

    Parameters
    ----------
    dataset: xr.Dataset
        The dataset containing the data variable
    data_var : xr.DataArray
        The data variable.
    freq : CLIMO_FREQ
        The frequency for calculating climatology.

    Returns
    -------
    xr.DataArray
        The variable's climatology.
    """
    ds = dataset.copy()
    ds = xc.center_times(ds)

    # The variable's time dim key is stored here for reuse in subsetting.
    time_dim = xc.get_dim_keys(ds[var_key], axis="T")

    if freq in ["ANN"]:
        ds_climo = ds.temporal.average(var_key, weighted=True)
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
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.month == int(freq))})
        ds_climo = ds.temporal.climatology(var_key, freq="month", weighted=True)
    elif freq in ["DJF", "MAM", "JJA", "SON"]:
        ds = ds.isel({f"{time_dim}": (ds[time_dim].dt.season == freq)})
        ds_climo = ds.temporal.climatology(var_key, freq="season", weighted=True)
    else:
        raise ValueError(
            f"`freq='{freq}'` is not a valid climatology frequency. Options "
            f"include {CLIMO_FREQS}'"
        )

    dv_climo = ds_climo[var_key].copy()

    # The time dimension should be a singleton after averaging so it should
    # be squeezed and dropped from the climatology data variable.
    if time_dim in dv_climo.dims:
        dv_climo = dv_climo.squeeze(dim=time_dim).drop_vars(time_dim)

    return dv_climo
