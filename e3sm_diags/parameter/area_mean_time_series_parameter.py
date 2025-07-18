from e3sm_diags.logger import _setup_child_logger

from .core_parameter import CoreParameter

logger = _setup_child_logger(__name__)


class AreaMeanTimeSeriesParameter(CoreParameter):
    def __init__(self):
        super(AreaMeanTimeSeriesParameter, self).__init__()
        # Override existing attributes
        # ============================
        self.ref_names = []
        self.ref_timeseries_input = True
        self.test_timeseries_input = True

        # Granulating with regions doesn't make sense, because we have multiple
        # regions for each plot. So keep all of the default values except
        # regions.
        self.granulate.remove("regions")
        self.granulate.remove("seasons")

        # Custom attributes
        # =================
        self.start_yr: str | None = None
        self.end_yr: str | None = None

    def check_values(self):
        if not self.ref_names:
            msg = "You have no value for ref_names. Calculate test data only"
            logger.info(msg)

        if self.start_yr is None and self.end_yr is None:
            msg = "You need to define both the 'start_yr' and 'end_yr' parameter."
            raise RuntimeError(msg)
