import xarray as xr

from e3sm_diags.logger import _setup_child_logger

logger = _setup_child_logger(__name__)


def monotonic(L):
    return _monotonically_increasing(L) or _monotonically_decreasing(L)


def _monotonically_decreasing(L):
    # FIXME: B905: zip() without an explicit strict= parameter
    return all(x >= y for x, y in zip(L, L[1:], strict=False))


def _monotonically_increasing(L):
    # FIXME: B905: zip() without an explicit strict= parameter
    return all(x <= y for x, y in zip(L, L[1:], strict=False))


def pad_year(year: int | str) -> str:
    """Pad the year with leading zeros to ensure it is 4 digits.

    This function ensures that the input year is properly formatted as a
    4-digit string, which is required for ISO-8601 date formats (YYYY-MM-DD).
    If the year is less than 1000, it is padded with leading zeros.

    Parameters
    ----------
    year : int or str
        The year to pad. Must be a non-negative integer or a string representing
        a non-negative integer. Floats are not allowed.

    Returns
    -------
    str
        The padded year as a 4-digit string (e.g., "0042" for year 42).

    Raises
    ------
    ValueError
        If the input year is not a non-negative integer or a string representing
        a non-negative integer, or if it is a float, or if it is outside the
        range 0 to 9999 inclusive.

    Examples
    --------
    >>> pad_year(42)
    '0042'
    >>> pad_year("42")
    '0042'
    >>> pad_year(2023)
    '2023'
    >>> pad_year("2023")
    '2023'
    >>> pad_year(-1)
    Traceback (most recent call last):
        ...
    ValueError: Year must be between 0 and 9999 inclusive.
    >>> pad_year(10000)
    Traceback (most recent call last):
        ...
    ValueError: Year must be between 0 and 9999 inclusive.
    >>> pad_year(42.0)
    Traceback (most recent call last):
        ...
    ValueError: Year must not be a float.
    """
    try:
        if isinstance(year, float):
            raise ValueError("Year must not be a float.")

        year = int(year)
    except (ValueError, TypeError) as e:
        raise ValueError(
            "Year must be a non-negative integer or a string representing a "
            "non-negative integer."
        ) from e

    if year < 0 or year > 9999:
        raise ValueError("Year must be between 0 and 9999 inclusive.")

    return f"{year:04d}"


def subtract_dataarrays(a: xr.DataArray, b: xr.DataArray) -> xr.DataArray:
    """Subtract two xarray DataArrays, preserving attributes from the left operand.

    This preserves the xarray <2025.11.0 behavior of preserving attributes from
    the left operand.

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
