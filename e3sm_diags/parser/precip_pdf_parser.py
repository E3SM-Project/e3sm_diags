from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.parser.core_parser import CoreParser


class PrecipPDFParser(CoreParser):
    def __init__(self, *args, **kwargs):
        kwargs["parameter_cls"] = PrecipPDFParameter
        super().__init__(*args, **kwargs)

    def add_arguments(self):
        super().add_arguments()

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
            "--regions",
            dest="regions",
            help="Regions for PDF calculation (e.g., 'tropics', 'conus').",
            nargs="+",
            required=False,
        )

        self.parser.add_argument(
            "--season_subset",
            dest="season_subset",
            help="Generate PDFs for all seasons (DJF, MAM, JJA, SON) in addition to annual. Default is False.",
            required=False,
        )
