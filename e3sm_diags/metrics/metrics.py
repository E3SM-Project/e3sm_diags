"""This module stores functions to calculate metrics using Xarray objects."""
from __future__ import annotations

from typing import List

import xarray as xr
import xcdat as xc
import xskillscore as xs

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)

AXES = ["X", "Y"]


def get_weights(ds: xr.Dataset):
    """Get weights for the X and Y spatial axes.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.

    Returns
    -------
    xr.DataArray
        Weights for the specified axis.
    """
    return ds.spatial.get_weights(axis=["X", "Y"])


def spatial_avg(
    ds: xr.Dataset, var_key: str, as_list: bool = True
) -> List[float] | xr.DataArray:
    """Compute a variable's weighted spatial average.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the varible.
    as_list : bool
        Return the spatial average as a list of floats, by default True.
        If False, return an xr.DataArray.

    Returns
    -------
    List[float] | xr.DataArray
        The spatial average of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.

    Notes
    -----
    Replaces `e3sm_diags.metrics.mean`.
    """
    ds_avg = ds.spatial.average(var_key, axis=AXES, weights="generate")
    results = ds_avg[var_key]

    if as_list:
        return results.data.tolist()

    return results


def std(ds: xr.Dataset, var_key: str) -> List[float]:
    """Compute the weighted standard deviation for a variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.

    Returns
    -------
    List[float]
        The standard deviation of the variable based on the specified axis.

    Raises
    ------
    ValueError
        If the axis argument contains an invalid value.

    Notes
    -----
    Replaces `e3sm_diags.metrics.std`.
    """
    dv = ds[var_key].copy()

    weights = ds.spatial.get_weights(axis=AXES, data_var=var_key)
    dims = _get_dims(dv, axis=AXES)

    result = dv.weighted(weights).std(dim=dims, keep_attrs=True)

    return result.data.tolist()


def correlation(ds_a: xr.Dataset, ds_b: xr.Dataset, var_key: str) -> List[float]:
    """Compute the correlation coefficient between two variables.

    This function uses the Pearson correlation coefficient. Refer to [1]_ for
    more information.

    Parameters
    ----------
    ds_a : xr.Dataset
        The first dataset.
    ds_b : xr.Dataset
        The second dataset.
    var_key: str
        The key of the variable.

    Returns
    -------
    List[float]
        The weighted correlation coefficient.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

    Notes
    -----
    Replaces `e3sm_diags.metrics.corr`.
    """
    var_a = ds_a[var_key]
    var_b = ds_b[var_key]

    # Dimensions, bounds, and coordinates should be identical between datasets,
    # so use the first dataset and variable to get dimensions and weights.
    dims = _get_dims(var_a, axis=AXES)
    weights = get_weights(ds_a)

    result = xs.pearson_r(var_a, var_b, dim=dims, weights=weights, skipna=True)
    results_list = result.data.tolist()

    return results_list


def rmse(ds_a: xr.Dataset, ds_b: xr.Dataset, var_key: str) -> List[float]:
    """Calculates the root mean square error (RMSE) between two variables.

    Parameters
    ----------
    ds_a : xr.Dataset
        The first dataset.
    ds_b : xr.Dataset
        The second dataset.
    var_key: str
        The key of the variable.

    Returns
    -------
    List[float]
        The root mean square error.

    Notes
    -----
    Replaces `e3sm_diags.metrics.rmse`.
    """
    var_a = ds_a[var_key]
    var_b = ds_b[var_key]

    # Dimensions, bounds, and coordinates should be identical between datasets,
    # so use the first dataset and variable to get dimensions and weights.
    dims = _get_dims(var_a, axis=AXES)
    weights = get_weights(ds_a)

    result = xs.rmse(var_a, var_b, dim=dims, weights=weights, skipna=True)
    results_list = result.data.tolist()

    return results_list


def _get_dims(da: xr.DataArray, axis: List[str]):
    """Get the dimensions for an axis in an xarray.DataArray.

    The dimensions are passed to the ``dim`` argument in xarray or xarray-based
    computational APIs, such as ``.std()``.

    Parameters
    ----------
    da : xr.DataArray
        The array.
    axis : List[str]
        A list of axis strings.

    Returns
    -------
    List[str]
        A list of dimensions.
    """
    dims = []

    for a in axis:
        dim_key = xc.get_dim_keys(da, axis=a)
        dims.append(dim_key)

    return dims
