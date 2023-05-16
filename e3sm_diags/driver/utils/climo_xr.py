"""The `climo` function from `climo.py` implemented to operate on xr.DataArray

NOTE: This function is being used as a temporary replacement of `climo.py`
for refactoring the codebase with xarray. Eventually, `climo_xcdat.py` will
be used instead (refer to that module's docstring for more info).
"""
from typing import Dict, List, Literal, get_args

import numpy as np
import numpy.ma as ma
import xarray as xr
import xcdat as xc

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)

# A type annotation and list of accepted climatology frequencies.
CLIMO_FREQ = Literal[
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
    "ANN",
    "DJF",
    "MAM",
    "JJA",
    "SON",
]
CLIMO_FREQS = get_args(CLIMO_FREQ)

# A dictionary mapping climatology frequencies to their months or seasons.
CLIMO_CYCLE_MAP = {
    "ANNUALCYCLE": [
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
    ],
    "SEASONALCYCLE": ["DJF", "MAM", "JJA", "SON"],
}
# A dictionary mapping climatology frequencies to their indexes for grouping
# coordinate points for weighted averaging.
FREQ_IDX_MAP: Dict[CLIMO_FREQ, List[int]] = {
    "01": [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "02": [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "03": [0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0],
    "04": [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
    "05": [0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0],
    "06": [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
    "07": [0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0],
    "08": [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
    "09": [0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0],
    "10": [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    "11": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0],
    "12": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "DJF": [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
    "MAM": [0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0],
    "JJA": [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0],
    "SON": [0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 0],
    "ANN": [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1],
}


def climo(data_var: xr.DataArray, freq: CLIMO_FREQ):
    """Computes a variable's climatology for the given frequency.

    The data variable's source dataset is retrieved for centering time
    coordinates using the time bounds. If there is no "source" filepath
    associated with the data variable, then it is considered a derived
    variable (a variable created from other variables in the dataset). In this
    case, the data variable is converted to a dataset to add time bounds and
    center the time coordinates.

    After averaging, the data variable's time dimension is squeezed and the
    time coordinates because they become singletons.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.
    freq : CLIMO_FREQ
        The frequency for calculating climatology.

    Returns
    -------
    xr.DataArray
        The variable's climatology.
    """
    ds = _get_ds_from_dv(data_var)

    try:
        dv_time = xc.get_dim_coords(data_var, axis="T")
    except KeyError:
        logger.warning(
            f"The climatology for {data_var.key} could not be calculated because it "
            "does not have time coordinates"
        )
        return data_var

    # Compute time lengths to use for weights.
    time_bnds = ds.bounds.get_bounds(axis="T")
    time_lengths = (time_bnds[:, 1] - time_bnds[:, 0]).astype(np.float64)

    # Get the frequency's cycle index map and number of cycles.
    if freq not in get_args(CLIMO_FREQ):
        raise ValueError(
            f"`freq='{freq}'` is not a valid climatology frequency. Options "
            f"include {get_args(CLIMO_FREQ)}'"
        )

    # Convert data variable from an `xr.DataArray` to a `np.MaskedArray` to
    # utilize the weighted averaging function.
    dv_masked = data_var.to_masked_array()

    # The indexes of the time coordinates related to the frequency are retrieved
    # using the frequency index map (`FREQ_IDX_MAP``). The data variable is
    # subsetted using the time indexes before calculating its weighted average.
    time_idx = []
    for i in range(len(dv_time)):
        month = dv_time[i].dt.month.item()
        idx = FREQ_IDX_MAP[freq][month - 1]
        time_idx.append(idx)

    time_idx = np.array(time_idx, dtype=np.int64).nonzero()
    climo = ma.average(dv_masked[time_idx], axis=0, weights=time_lengths[time_idx])

    # Construct the climatology xr.DataArray using the climatology output.
    # The time dimension/coords are not included since it is squeezed during
    # averaging.
    dims = [dim for dim in data_var.dims if dim != dv_time.name]
    coords = {k: v for k, v in data_var.coords.items() if k in dims}
    dv_climo = xr.DataArray(
        name=data_var.name,
        data=climo,
        coords={**coords},
        dims=dims,
        attrs=data_var.attrs,
    )

    return dv_climo


def _get_ds_from_dv(data_var: xr.DataArray) -> xr.Dataset:
    """Get the data variable's source Dataset.

    If there is no "source" filepath associated with the data variable, then it
    is considered a derived variable (a variable created from other variables
    in the dataset). We need to convert this data variable to a dataset to
    center the time coordinates and add time bounds.

    Parameters
    ----------
    data_var : xr.DataArray
        The data variable.

    Returns
    -------
    xr.Dataset
        The source dataset.
    """
    filepath = data_var.encoding.get("source")

    if filepath is None:
        ds = data_var.to_dataset()
        ds = xc.center_times(ds)
        ds = ds.bounds.add_bounds(axis="T")
    else:
        ds = xc.open_dataset(filepath, center_times=True)

    return ds
