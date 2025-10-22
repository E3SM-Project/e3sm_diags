from __future__ import annotations

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

        # Style options for native grid visualization
        self.edge_color = None  # Set to a color string to show grid edges
        self.edge_width = 0.3  # Width of grid edges when displayed

        # Option to disable the grid antialiasing (may improve performance)
        self.antialiased = False

    def check_values(self):
        """Verifies that required values are properly set."""
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
            super()._set_name_yrs_attrs(test_ds, ref_ds, season)
        else:
            self._set_time_slice_name_yrs_attrs(test_ds, ref_ds, season)

    def _set_time_slice_name_yrs_attrs(
        self, test_ds: Dataset, ref_ds: Dataset, time_slice: TimeSlice
    ) -> None:
        """Set name_yrs attributes for time slice-based processing.

        This function sets up the necessary attributes for file naming and
        processing when using time_slices instead of seasons.

        NOTE: This method is only used by lat_lon_native. Other diagnostic sets
        call self._set_name_yrs_attrs() from the parent class directly, even
        with time slices. We may want to refactor this in the future for
        consistency.

        Parameters
        ----------
        test_ds : Dataset
            The test dataset object.
        ref_ds : Dataset
            The reference dataset object.
        time_slice : TimeSlice
            The time slice specification.

        Notes
        -----
        This function modifies the parameter object in-place by setting:
        - parameter.current_time_slice
        - parameter.test_name_yrs
        - parameter.ref_name_yrs
        """
        # Set the time slice info for potential use in plotting/output
        self.current_time_slice = time_slice

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
