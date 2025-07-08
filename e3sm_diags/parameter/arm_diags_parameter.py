from .core_parameter import CoreParameter


class ARMDiagsParameter(CoreParameter):
    def __init__(self):
        super(ARMDiagsParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.granulate.remove("seasons")
        self.test_timeseries_input = True
        self.ref_timeseries_input = True

        # Custom attributes
        # =============================
        self.test_start_yr: str | None = None
        self.test_end_yr: str | None = None
        self.ref_start_yr: str | None = None
        self.ref_end_yr: str | None = None

        # "Options include: annual_cycle", "diurnal_cycle", "diurnal_cycle_zt",
        # "pdf_daily", "convection_onset"
        self.diags_set: str | None = None

        self.var_name: str | None = None
        self.var_units: str | None = None

        # The time interval in hours for the diurnal cycle
        self.time_interval: int | None = None
