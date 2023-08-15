from __future__ import annotations

from typing import List

import xarray as xr
import xcdat as xc
import xskillscore as xs

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def get_weights(ds: xr.Dataset, axis: List[str] = ["X", "Y"]):
    """Get axis weights for an axis/axes in the dataset.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset.
    axis : List[str], optional
        A list of axis, by default ["X", "Y"].

    Returns
    -------
    xr.DataArray
        Weights for the specified axis.
    """
    _validate_axis_arg(axis)

    return ds.spatial.get_weights(axis=axis)


def spatial_avg(
    ds: xr.Dataset, var_key: str, axis: List[str] = ["X", "Y"], serialize: bool = False
) -> xr.DataArray | List[float]:
    """Compute a variable's weighted spatial average.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the varible.
    axis : List[str], optional
        The axis to compute spatial average on, by default ["X", "Y"]. Options
        include "X" and "Y".
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
    This function is intended to replace ``e3sm_diags.metrics.mean()``.
    """
    _validate_axis_arg(axis)

    ds_avg = ds.spatial.average(var_key, axis=axis, weights="generate")
    results = ds_avg[var_key]

    if serialize:
        return results.data.tolist()

    return results


def std(
    ds: xr.Dataset, var_key: str, axis=["X", "Y"], serialize: bool = False
) -> xr.DataArray | List[float]:
    """Compute the weighted standard deviation for a variable.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset containing the variable.
    var_key : str
        The key of the variable.
    axis : List[str], optional
        The spatial axis to compute standard deviation on, by default
        ["X", "Y"]. Options include "X" and "Y".
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
    This function is intended to replace ``e3sm_diags.metrics.std()``.
    """
    _validate_axis_arg(axis)

    dv = ds[var_key].copy()

    weights = ds.spatial.get_weights(axis, data_var=var_key)
    dims = _get_dims(dv, axis)

    result = dv.weighted(weights).std(dim=dims, keep_attrs=True)

    if result:
        return result.data.tolist()

    return result


def correlation(
    da_a: xr.DataArray,
    da_b: xr.DataArray,
    axis: List[str] = ["X", "Y"],
    weights: xr.DataArray = None,
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
    axis : List[str] , optional
        The axis to compute the correlation on, by default ["X", "Y"]
    weights: xr.DataArray, optional
        The weight related to the specified ``axis``, by default None.
        If None, the results are unweighted.
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
    This function is intended to replace ``e3sm_diags.metrics.corr()``.
    """
    _validate_axis_arg(axis)
    dims = _get_dims(da_a, axis)

    result = xs.pearson_r(da_a, da_b, dim=dims, weights=weights)

    if serialize:
        return result.data.tolist()

    return result


def rmse(
    da_a: xr.DataArray,
    da_b: xr.DataArray,
    axis: List[str] = ["X", "Y"],
    weights: xr.DataArray = None,
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
    weights: xr.DataArray, optional
        The weight related to the specified ``axis``, by default None.
        If None, the results are unweighted.
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
    This function is intended to replace ``e3sm_diags.metrics.rmse()``.
    """
    _validate_axis_arg(axis)
    dims = _get_dims(da_a, axis)

    result = xs.rmse(da_a, da_b, dim=dims, weights=weights)

    if serialize:
        return result.data.tolist()

    return result


def _validate_axis_arg(axis: List[str]):
    """Validates the ``axis`` argument is a list with supported values.

    Supported values include "X", "Y", and "T".

    Parameters
    ----------
    axis : List[str]
        A list of axis strings.

    Raises
    ------
    ValueError
        If ``axis`` contains an unsupported value(s).
    """
    for k in axis:
        if k not in ["X", "Y"]:
            raise ValueError(
                f"The `axis` argument has an unsupported value ('{k}'). "
                "Supported values include: ['X'], ['Y'], ['X', 'Y']."
            )


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
