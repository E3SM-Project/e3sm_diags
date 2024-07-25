"""This module stores climatology functions operating on Xarray objects.

NOTE: Replaces `e3sm_diags.driver.utils.climo`.

This file will eventually be refactored to use xCDAT's climatology API.
"""
from typing import Dict, List, Literal, get_args

import numpy as np
import numpy.ma as ma
import xarray as xr
import xcdat as xc

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)

# A type annotation and list representing accepted climatology frequencies.
# Accepted frequencies include the month integer and season string.
ClimoFreq = Literal[
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
    "ANNUALCYCLE",
]
CLIMO_FREQS = get_args(ClimoFreq)

# A dictionary that maps climatology frequencies to the appropriate cycle
# for grouping.
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
FREQ_IDX_MAP: Dict[ClimoFreq, List[int]] = {
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


def climo(dataset: xr.Dataset, var_key: str, freq: ClimoFreq):
    """Computes a variable's climatology for the given frequency.

    Parameters
    ----------
    dataset: xr.Dataset
        The time series dataset.
    var_key : xr.DataArray
        The key of the variable in the Dataset to calculate climatology for.
    freq : CLIMO_FREQ
        The frequency for calculating climatology.

    Returns
    -------
    xr.DataArray
        The variable's climatology.
    """
    # Get the frequency's cycle index map and number of cycles.
    if freq not in get_args(ClimoFreq):
        raise ValueError(
            f"`freq='{freq}'` is not a valid climatology frequency. Options "
            f"include {get_args(ClimoFreq)}'"
        )

    # Time coordinates are centered (if they aren't already) for more robust
    # weighted averaging calculations.
    ds = dataset.copy()
    ds = xc.center_times(ds)

    # Extract the data variable from the new dataset to calculate weighted
    # averaging.
    dv = ds[var_key].copy()
    time_coords = xc.get_dim_coords(dv, axis="T")

    # Loop over the time coordinates to get the indexes related to the
    # user-specified climatology frequency using the frequency index map
    # (`FREQ_IDX_MAP``).
    time_idx = []
    for i in range(len(time_coords)):
        month = time_coords[i].dt.month.item()
        idx = FREQ_IDX_MAP[freq][month - 1]
        time_idx.append(idx)

    time_idx = np.array(time_idx, dtype=np.int64).nonzero()  # type: ignore

    # Convert data variable from an `xr.DataArray` to a `np.MaskedArray` to
    # utilize the weighted averaging function and use the time bounds
    # to calculate time lengths for weights.
    # NOTE: Since `time_bnds`` are decoded, the arithmetic to produce
    # `time_lengths` will result in the weighted averaging having an extremely
    # small floating point difference (1e-16+) compared to `climo.py`.
    dv_masked = dv.to_masked_array()

    time_bnds = ds.bounds.get_bounds(axis="T")
    time_lengths = (time_bnds[:, 1] - time_bnds[:, 0]).astype(np.float64)

    # Calculate the weighted average of the masked data variable using the
    # appropriate indexes and weights.
    climo = ma.average(dv_masked[time_idx], axis=0, weights=time_lengths[time_idx])

    # Construct the climatology xr.DataArray using the averaging output. The
    # time coordinates are not included since they become a singleton after
    # averaging.
    dims = [dim for dim in dv.dims if dim != time_coords.name]
    coords = {k: v for k, v in dv.coords.items() if k in dims}
    dv_climo = xr.DataArray(
        name=dv.name,
        data=climo,
        coords={**coords},
        dims=dims,
        attrs=dv.attrs,
    )

    return dv_climo
