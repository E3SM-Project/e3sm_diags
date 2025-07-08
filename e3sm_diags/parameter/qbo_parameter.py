from .time_series_parameter import TimeSeriesParameter


class QboParameter(TimeSeriesParameter):
    def __init__(self):
        super(QboParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.print_statements = False
        self.ref_timeseries_input = True
        self.test_timeseries_input = True
        self.granulate.remove("seasons")

        # Custom attributes
        # -----------------
        self.ref_yrs: str | None = None
        self.test_yrs: str | None = None
