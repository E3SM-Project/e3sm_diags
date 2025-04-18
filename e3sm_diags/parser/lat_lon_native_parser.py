from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter
from e3sm_diags.parser.core_parser import CoreParser


class LatLonNativeParser(CoreParser):
    def __init__(self, *args, **kwargs):
        # FIXME: B026 Star-arg unpacking after a keyword argument is strongly discouraged
        super().__init__(parameter_cls=LatLonNativeParameter, *args, **kwargs)  # type: ignore # noqa: B026

    def add_arguments(self):
        super().add_arguments()

        self.parser.add_argument(
            "--grid_file",
            dest="grid_file",
            help="Path to the native grid file to use for visualization",
            required=False,
        )

        self.parser.add_argument(
            "--split_periodic_elements",
            dest="split_periodic_elements",
            help="Whether to split elements that cross the date line for better visualization",
            action="store_true",
            default=True,
            required=False,
        )

        self.parser.add_argument(
            "--antialiased",
            dest="antialiased",
            help="Whether to use antialiasing for grid edges",
            action="store_true",
            default=False,
            required=False,
        )

        self.parser.add_argument(
            "--edge_color",
            dest="edge_color",
            help="Color for grid edges (None for no edges)",
            required=False,
        )

        self.parser.add_argument(
            "--edge_width",
            dest="edge_width",
            type=float,
            default=0.3,
            help="Width of grid edges when displayed",
            required=False,
        )
