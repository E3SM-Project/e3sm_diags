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
    ds: xr.Dataset, var_key: str, serialize: bool = False
) -> xr.DataArray | List[float]:
    """Compute a variable's weighted spatial average.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the varible.
    serialize : bool, optional
        If True, convert the underlying `np.array` to a Python list, by default
        False. `np.float32` and `np.float64` arrays are not JSON serializable
        and must be converted to a Python list of native floats if dumping
        the function's output to a JSON file.

    Returns
    -------
    xr.DataArray | List[float]
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

    if serialize:
        return results.data.tolist()

    return results


def std(
    ds: xr.Dataset, var_key: str, serialize: bool = False
) -> xr.DataArray | List[float]:
    """Compute the weighted standard deviation for a variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    serialize : bool, optional
        If True, convert the underlying `np.array` to a Python list, by default
        False. `np.float32` and `np.float64` arrays are not JSON serializable
        and must be converted to a Python list of native floats if dumping
        the function's output to a JSON file.

    Returns
    -------
    xr.DataArray | List[float]
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

    if serialize:
        return result.data.tolist()

    return result


def correlation(
    da_a: xr.DataArray,
    da_b: xr.DataArray,
    weights: xr.DataArray,
    serialize: bool = False,
) -> xr.DataArray:
    """Compute the correlation coefficient between two variables.

    This function uses the Pearson correlation coefficient. Refer to [1]_ for
    more information.

    Parameters
    ----------
    da_a : xr.DataArray
        The first variable.
    da_b : xr.DataArray
        The second variable.
    weights: xr.DataArray
        The weights for the X and Y axes.
    serialize : bool, optional
        If True, convert the underlying `np.array` to a Python list, by default
        False. `np.float32` and `np.float64` arrays are not JSON serializable
        and must be converted to a Python list of native floats if dumping
        the function's output to a JSON file.

    Returns
    -------
    xr.DataArray
        The weighted correlation coefficient.

    References
    ----------

    .. [1] https://en.wikipedia.org/wiki/Pearson_correlation_coefficient

    Notes
    -----
    Replaces `e3sm_diags.metrics.corr`.
    """
    dims = _get_dims(da_a, axis=["X", "Y"])

    result = xs.pearson_r(da_a, da_b, dim=dims, weights=weights, skipna=True)

    if serialize:
        return result.data.tolist()

    return result


def rmse(
    da_a: xr.DataArray,
    da_b: xr.DataArray,
    weights: xr.DataArray,
    serialize: bool = False,
) -> xr.DataArray | List[float]:
    """Calculates the root mean square error (RMSE) between two variables.

    Parameters
    ----------
    da_a : xr.DataArray
        The first variable.
    da_b : xr.DataArray
        The second variable.
    axis : List[str] , optional
        The axis to compute the correlation on, by default ["X", "Y"]
    weights: xr.DataArray
        The weights for the X and Y axes.
    serialize : bool, optional
        If True, convert the underlying `np.array` to a Python list, by default
        False. `np.float32` and `np.float64` arrays are not JSON serializable
        and must be converted to a Python list of native floats if dumping
        the function's output to a JSON file.

    Returns
    -------
    xr.DataArray | List[float]
        The root mean square error.

    Notes
    -----
    Replaces `e3sm_diags.metrics.rmse`.
    """
    dims = _get_dims(da_a, axis=AXES)

    result = xs.rmse(da_a, da_b, dim=dims, weights=weights, skipna=True)

    if serialize:
        return result.data.tolist()

    return result


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
