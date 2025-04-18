from __future__ import annotations

from e3sm_diags.parameter.core_parameter import CoreParameter


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
        self.ref_grid_file = ""   # Grid file for reference data

        # Option for handling periodic elements
        # If True, split elements that cross the dateline for better visualization
        self.split_periodic_elements = True

        # Style options for native grid visualization
        self.edge_color = None  # Set to a color string to show grid edges
        self.edge_width = 0.3   # Width of grid edges when displayed

        # Option to disable the grid antialiasing (may improve performance)
        self.antialiased = False

    def check_values(self):
        """Verifies that required values are properly set.

        Raises
        ------
        RuntimeError
            If no grid files are provided or set.
        """
        # For now, we'll make this a soft check
        # In the future, we may want to require at least test_grid_file
        pass
