from typing import Optional

from .time_series_parameter import TimeSeriesParameter


class PrecipPDFParameter(TimeSeriesParameter):
    def __init__(self):
        super(PrecipPDFParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        self.granulate.remove("seasons")

        # Custom attributes
        # -----------------
        self.ref_yrs: Optional[str] = None
        self.test_yrs: Optional[str] = None

        # Regions for PDF calculation (use built-in region names from REGION_SPECS)
        # Default regions match tropical subseasonal analysis
        self.regions = ["15S15N"]  # Default: tropical region (15S-15N)

        # PDF bins (consistent with GPCP/TRMM)
        self.num_bins = 129  # Number of precipitation bins

        # Season subsetting (different from seasonal averaging in granulate)
        # When True: generate PDFs for all months (annual) plus each season (DJF, MAM, JJA, SON)
        # When False: only generate annual PDF (all months)
        self.season_subset: bool = True
