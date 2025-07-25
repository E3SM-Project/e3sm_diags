from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parser.core_parser import CoreParser


class EnsoDiagsParser(CoreParser):
    def __init__(self, *args, **kwargs):
        # FIXME: B026 Star-arg unpacking after a keyword argument is strongly discouraged
        super().__init__(parameter_cls=EnsoDiagsParameter, *args, **kwargs)  # type: ignore  # noqa: B026

    def add_arguments(self):
        super().add_arguments()

        self.parser.add_argument(
            "--ref_names",
            type=str,
            nargs="+",
            dest="ref_names",
            help="List of reference names.",
            required=False,
        )

        self.parser.add_argument(
            "--ref_timeseries_input",
            dest="ref_timeseries_input",
            help="The input reference data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--test_timeseries_input",
            dest="test_timeseries_input",
            help="The input test data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.parser.add_argument(
            "--start_yr",
            dest="start_yr",
            help="Start year for the timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--end_yr",
            dest="end_yr",
            help="End year for the timeseries files.",
            required=False,
        )

        self.parser.add_argument(
            "--plot_type",
            dest="plot_type",
            help="Type of plot to generate: 'map' or 'scatter'.",
            required=False,
        )
