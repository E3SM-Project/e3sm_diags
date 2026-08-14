from .time_series_parameter import TimeSeriesParameter


class EnsoDiagsParameter(TimeSeriesParameter):
    def __init__(self):
        super(EnsoDiagsParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.granulate.remove("seasons")
        self.ref_timeseries_input = True
        self.test_timeseries_input = True

        # Custom attributes
        # =============================
        self.nino_region = "NINO34"
        # The list of nino regions stacked (one subplot row each) in the
        # "nino_index_timeseries" and "seasonality" plot types.
        self.nino_regions = ["NINO3", "NINO34", "NINO4"]
        # The lags (in months) tiled as subplot rows in the "lead_lag" plot
        # type. Positive lags indicate the nino index leading the field.
        self.lead_lag_months = [-8, -4, 0, 4, 8]
        # Figures produced by the "lead_lag" plot type, recorded by the driver
        # so the viewer can link each one (a single run makes both the
        # regression and correlation figures). Each entry is a dict with
        # "output_file" and "descr" keys.
        self.lead_lag_entries: list[dict] = []
        self.plot_type = "regression_map"
        self.print_statements = False

    def check_values(self):
        super(EnsoDiagsParameter, self).check_values()
        valid_nino_regions = ["NINO3", "NINO34", "NINO4"]
        if self.nino_region not in valid_nino_regions:
            msg = "nino_region={} not in {}".format(
                self.nino_region, valid_nino_regions
            )
            raise RuntimeError(msg)

        valid_plot_types = [
            "regression_map",
            "feedback",
            "nino_index_timeseries",
            "seasonality",
            "interannual_variability",
            "equatorial_soi",
            "lead_lag",
        ]
        if self.plot_type not in valid_plot_types:
            msg = "plot_type={} not in {}".format(self.plot_type, valid_plot_types)
            raise RuntimeError(msg)
