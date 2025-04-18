from __future__ import annotations

from e3sm_diags.parameter.core_parameter import CoreParameter


class LatLonNativeParameter(CoreParameter):
    """
    Parameters for the lat_lon_native diagnostic set.
    """

    def __init__(self):
        super().__init__()

        # Path to the grid file for the native grid
        self.grid_file = ""

        # This field is deprecated - use grid_file directly

        # Option for handling periodic elements
        # If True, split elements that cross the dateline for better visualization
        self.split_periodic_elements = True

        # Style options for native grid visualization
        self.edge_color = None  # Set to a color string to show grid edges
        self.edge_width = 0.3  # Width of grid edges when displayed

        # Option to disable the grid antialiasing (may improve performance)
        self.antialiased = False
