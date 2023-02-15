from e3sm_diags.parameter.qbo_parameter import QboParameter

from .core_parser_new import CoreParser


class QboParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=QboParameter, *args, **kwargs)  # type:ignore

    def load_default_args(self):

        self.add_argument(
            "--ref_timeseries_input",
            dest="ref_timeseries_input",
            help="The input reference data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--test_timeseries_input",
            dest="test_timeseries_input",
            help="The input test data are timeseries files.",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--start_yr",
            dest="start_yr",
            help="Start year for the timeseries files.",
            required=False,
        )

        self.add_argument(
            "--end_yr",
            dest="end_yr",
            help="End year for the timeseries files.",
            required=False,
        )
