from __future__ import annotations

import re
from typing import TYPE_CHECKING

from e3sm_diags.parameter.core_parameter import CoreParameter

if TYPE_CHECKING:
    from e3sm_diags.driver.utils.dataset_xr import Dataset
    from e3sm_diags.driver.utils.type_annotations import TimeSelection, TimeSlice


class LatLonNativeParameter(CoreParameter):
    """Parameters for the lat_lon_native diagnostic set.

    This diagnostic set allows displaying data on native grids (e.g., cubed-sphere)
    using uxarray's visualization capabilities.
    """

    def __init__(self):
        super(LatLonNativeParameter, self).__init__()

        # Override existing attributes
        # =============================
        # Path to the grid files for the native grids
        self.test_grid_file = ""  # Grid file for test data
        self.ref_grid_file = ""  # Grid file for reference data

        # Option for handling periodic elements
        # If True, split elements that cross the dateline for better visualization
        self.split_periodic_elements = True

        # Style options for native grid visualization
        self.edge_color = None  # Set to a color string to show grid edges
        self.edge_width = 0.3  # Width of grid edges when displayed

        # Option to disable the grid antialiasing (may improve performance)
        self.antialiased = False

        # Time selection parameters (mutually exclusive with seasons)
        # Either use seasons (inherited from CoreParameter) OR time_slices
        # Index-based time selection with stride support
        # Examples: ["0:10:2", "5:15", "7"] for start:end:stride, start:end, or single index
        self.time_slices: list[TimeSlice] = []

    def check_values(self):
        """Verifies that required values are properly set.

        Raises
        ------
        RuntimeError
            If no grid files are provided or set.
        RuntimeError
            If neither seasons nor time_slices are specified.
        """
        has_seasons = len(self.seasons) > 0
        has_time_slices = len(self.time_slices) > 0

        if not has_seasons and not has_time_slices:
            raise RuntimeError(
                "Must specify either 'seasons' or 'time_slices'. "
                "Use 'seasons' for climatological analysis (e.g., ['ANN', 'DJF']) "
                "or 'time_slices' for index-based selection (e.g., ['0:10:2', '5:15'])."
            )

        # Validate time_slice format if provided
        if has_time_slices:
            for time_slice in self.time_slices:
                self._validate_time_slice_format(time_slice)

        # TODO: For now, we'll make grid file check a soft check. In the future,
        # we may want to require at least test_grid_file
        pass

    def _validate_time_slice_format(self, time_slice: str):
        r"""Validate that time_slice follows the expected format.

        This regex pattern for slice notation is designed to match a
        latitude/longitude-like format with optional degrees, minutes, and
        seconds.
            - ^: Matches the start of the string.
            - (-?\d+|): Matches an optional integer (can be negative) for degrees.
            - (?::(-?\d+|): Matches an optional colon followed by an optional
            integer (can be negative) for minutes.
            - (?::(-?\d+|)): Matches an optional colon followed by an optional
            integer (can be negative) for seconds.
            - )?: Makes the minutes and seconds groups optional.
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
        """
        pattern = r"^(-?\d+|)(?::(-?\d+|)(?::(-?\d+|))?)?$"

        if not re.match(pattern, time_slice.strip()):
            raise ValueError(
                f"Invalid time_slice format: '{time_slice}'. "
                f"Expected formats: 'index', 'start:end', 'start:end:stride', "
                f"':end', 'start:', or '::stride'. Examples: '5', '0:10', '0:10:2'"
            )

    def _set_name_yrs_attrs(
        self, test_ds: Dataset, ref_ds: Dataset, season: TimeSelection | None
    ):
        """Override parent method to handle both ClimoFreq and time slice strings.

        Parameters
        ----------
        test_ds : Dataset
            The test dataset object.
        ref_ds : Dataset
            The reference dataset object.
        season : TimeSelection | None
            The season or time slice string.
        """
        from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQS

        if season is None or season in CLIMO_FREQS:
            # Standard climatology season, use parent implementation.
            super()._set_name_yrs_attrs(test_ds, ref_ds, season)
        else:
            # This is a time slice string, handle it specially.
            self._set_time_slice_attrs(test_ds, ref_ds, season)

    def _set_time_slice_attrs(self, test_ds: Dataset, ref_ds: Dataset, time_slice: str):
        """Set attributes for time slice-based processing.

        This method sets up the necessary attributes for file naming and
        processing when using time_slices instead of seasons.

        Store the time slice info but keep current_set as the diagnostic set name
        current_set should remain as "lat_lon_native" for proper directory structure
        The time slice will be used in filename generation via other attributes

        Parameters
        ----------
        test_ds : Dataset
            The test dataset object.
        ref_ds : Dataset
            The reference dataset object.
        time_slice : str
            The time slice specification.
        """
        # Set the time slice info for potential use in plotting/output
        self.current_time_slice = time_slice

        # For time slices, we manually set the name_yrs attributes instead of
        # calling parent method to avoid issues with the dataset's get_name_yrs_attr
        # expecting a valid season

        # Set test_name_yrs - use test dataset years if available, otherwise use
        # time slice info
        try:
            # Try to get year range from test dataset start/end years
            if hasattr(test_ds, "start_yr") and hasattr(test_ds, "end_yr"):
                test_years = f"{test_ds.start_yr:04d}-{test_ds.end_yr:04d}"
            else:
                test_years = f"timeslice_{time_slice}"

            self.test_name_yrs = f"{getattr(self, 'test_name', 'test')}_{test_years}"
        except AttributeError:
            self.test_name_yrs = (
                f"{getattr(self, 'test_name', 'test')}_timeslice_{time_slice}"
            )

        # Set ref_name_yrs - use ref dataset years if available, otherwise use time slice info
        try:
            if hasattr(ref_ds, "start_yr") and hasattr(ref_ds, "end_yr"):
                ref_years = f"{ref_ds.start_yr:04d}-{ref_ds.end_yr:04d}"
            else:
                ref_years = f"timeslice_{time_slice}"

            self.ref_name_yrs = f"{getattr(self, 'ref_name', 'ref')}_{ref_years}"
        except AttributeError:
            self.ref_name_yrs = (
                f"{getattr(self, 'ref_name', 'ref')}_timeslice_{time_slice}"
            )
