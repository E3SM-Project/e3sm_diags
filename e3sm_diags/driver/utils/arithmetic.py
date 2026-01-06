import xarray as xr


def subtract_dataarrays(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
    """Subtract two xarray DataArrays, preserving attributes from the left operand.

    This replicates the behavior of xarray <2025.11.0 where attributes are
    preserved from the left operand.

    Parameters
    ----------
    a : xr.DataArray
        The first DataArray.
    b : xr.DataArray
        The second DataArray.

    Returns
    -------
    xr.DataArray
        The result of a - b.
    """
    if a.shape != b.shape:
        raise ValueError("Input DataArrays must have the same shape.")

    result = a - b
    result.attrs = a.attrs.copy()

    return result


def add_dataarrays(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
    """Add two xarray DataArrays, preserving attributes from the left operand.

    This replicates the behavior of xarray <2025.11.0 where attributes are
    preserved from the left operand.

    Parameters
    ----------
    a : xr.DataArray
        The first DataArray.
    b : xr.DataArray
        The second DataArray.

    Returns
    -------
    xr.DataArray
        The result of a + b.
    """
    if a.shape != b.shape:
        raise ValueError("Input DataArrays must have the same shape.")

    result = a + b
    result.attrs = a.attrs.copy()

    return result
