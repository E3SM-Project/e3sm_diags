from __future__ import annotations

from typing import TYPE_CHECKING

from e3sm_diags.driver.utils.time_slice import (
    check_time_selection,
    set_time_slice_name_yrs_attrs,
    validate_time_slice_format,
)
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

        # Style options for native grid visualization
        self.edge_color = None  # Set to a color string to show grid edges
        self.edge_width = 0.3  # Width of grid edges when displayed

        # Option to disable the grid antialiasing (may improve performance)
        self.antialiased = False

        # Time selection parameters (mutually exclusive with seasons)
        # Either use seasons (inherited from CoreParameter) OR time_slices
        # Index-based time selection for snapshot analysis using individual time indices
        # Examples: ["0"], ["5"], ["0", "1", "2"]
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
        has_seasons, has_time_slices = check_time_selection(
            self.seasons, self.time_slices, require_one=True
        )

        # Validate time_slice format if provided
        if has_time_slices:
            for time_slice in self.time_slices:
                validate_time_slice_format(time_slice)

        # TODO: For now, we'll make grid file check a soft check. In the future,
        # we may want to require at least test_grid_file
        pass

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
        # Use shared utility function to set time slice attributes
        set_time_slice_name_yrs_attrs(self, test_ds, ref_ds, time_slice)
