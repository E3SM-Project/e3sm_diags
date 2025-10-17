"""Shared utilities for time slice functionality across diagnostic sets.

This module provides common validation, attribute setting, and data loading
utilities for time slice support. It extracts shared logic from lat_lon_native
to enable consistent time slice behavior across multiple diagnostic sets.
"""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from e3sm_diags.driver.utils.dataset_xr import Dataset
    from e3sm_diags.driver.utils.type_annotations import TimeSlice


def validate_time_slice_format(time_slice: str) -> None:
    r"""Validate that time_slice follows the expected format.

    This regex pattern for slice notation is designed to match Python slice
    syntax with optional components:
        - ^: Matches the start of the string.
        - (-?\d+|): Matches an optional integer (can be negative) for start.
        - (?::(-?\d+|): Matches an optional colon followed by an optional
          integer (can be negative) for end.
        - (?::(-?\d+|)): Matches an optional colon followed by an optional
          integer (can be negative) for stride.
        - )?: Makes the end and stride groups optional.
        - $: Matches the end of the string.

    Valid formats:
        - "index" (single index): "5"
        - "start:end" (range): "0:10"
        - "start:end:stride" (range with stride): "0:10:2"
        - ":end" (from beginning): ":10"
        - "start:" (to end): "5:"
        - "::stride" (full range with stride): "::2"

    Parameters
    ----------
    time_slice : str
        The time slice string to validate

    Raises
    ------
    ValueError
        If the time slice format is invalid

    Examples
    --------
    >>> validate_time_slice_format("5")  # Single index
    >>> validate_time_slice_format("0:10")  # Range
    >>> validate_time_slice_format("0:10:2")  # Range with stride
    >>> validate_time_slice_format(":10")  # From beginning
    >>> validate_time_slice_format("5:")  # To end
    >>> validate_time_slice_format("::2")  # Full range with stride
    """
    pattern = r"^(-?\d+|)(?::(-?\d+|)(?::(-?\d+|))?)?$"

    if not re.match(pattern, time_slice.strip()):
        raise ValueError(
            f"Invalid time_slice format: '{time_slice}'. "
            f"Expected formats: 'index', 'start:end', 'start:end:stride', "
            f"':end', 'start:', or '::stride'. Examples: '5', '0:10', '0:10:2'"
        )


def set_time_slice_name_yrs_attrs(
    parameter,
    test_ds: Dataset,
    ref_ds: Dataset,
    time_slice: str,
) -> None:
    """Set name_yrs attributes for time slice-based processing.

    This function sets up the necessary attributes for file naming and
    processing when using time_slices instead of seasons.

    Parameters
    ----------
    parameter : CoreParameter or subclass
        The parameter object to update.
    test_ds : Dataset
        The test dataset object.
    ref_ds : Dataset
        The reference dataset object.
    time_slice : str
        The time slice specification.

    Notes
    -----
    This function modifies the parameter object in-place by setting:
    - parameter.current_time_slice
    - parameter.test_name_yrs
    - parameter.ref_name_yrs
    """
    # Set the time slice info for potential use in plotting/output
    parameter.current_time_slice = time_slice

    # Set test_name_yrs - use test dataset years if available, otherwise use
    # time slice info
    try:
        # Try to get year range from test dataset start/end years
        if hasattr(test_ds, "start_yr") and hasattr(test_ds, "end_yr"):
            test_years = f"{test_ds.start_yr:04d}-{test_ds.end_yr:04d}"
        else:
            test_years = f"timeslice_{time_slice}"

        parameter.test_name_yrs = (
            f"{getattr(parameter, 'test_name', 'test')}_{test_years}"
        )
    except AttributeError:
        parameter.test_name_yrs = (
            f"{getattr(parameter, 'test_name', 'test')}_timeslice_{time_slice}"
        )

    # Set ref_name_yrs - use ref dataset years if available, otherwise use time slice info
    try:
        if hasattr(ref_ds, "start_yr") and hasattr(ref_ds, "end_yr"):
            ref_years = f"{ref_ds.start_yr:04d}-{ref_ds.end_yr:04d}"
        else:
            ref_years = f"timeslice_{time_slice}"

        parameter.ref_name_yrs = f"{getattr(parameter, 'ref_name', 'ref')}_{ref_years}"
    except AttributeError:
        parameter.ref_name_yrs = (
            f"{getattr(parameter, 'ref_name', 'ref')}_timeslice_{time_slice}"
        )


def check_time_selection(
    seasons: list,
    time_slices: list[TimeSlice],
    require_one: bool = True,
) -> tuple[bool, bool]:
    """Check time selection parameters (seasons vs time_slices).

    Parameters
    ----------
    seasons : list
        List of season specifications.
    time_slices : list[TimeSlice]
        List of time slice specifications.
    require_one : bool, optional
        If True, requires at least one of seasons or time_slices to be specified.
        Default is True.

    Returns
    -------
    tuple[bool, bool]
        A tuple of (has_seasons, has_time_slices)

    Raises
    ------
    RuntimeError
        If require_one is True and neither seasons nor time_slices are specified.
    """
    has_seasons = len(seasons) > 0
    has_time_slices = len(time_slices) > 0

    if require_one and not has_seasons and not has_time_slices:
        raise RuntimeError(
            "Must specify either 'seasons' or 'time_slices'. "
            "Use 'seasons' for climatological analysis (e.g., ['ANN', 'DJF']) "
            "or 'time_slices' for index-based selection (e.g., ['0:10:2', '5:15'])."
        )

    return has_seasons, has_time_slices
