import numpy

from e3sm_diags.driver.utils.general import monotonic

from .core_parameter import CoreParameter


class ZonalMean2dStratosphereParameter(CoreParameter):
    def __init__(self):
        super(ZonalMean2dStratosphereParameter, self).__init__()
        self.plevs = numpy.logspace(0, 2.0, num=10).tolist()
        self.plot_log_plevs = True
        self.plot_plevs = False
        # Granulating plevs causes duplicate plots in this case.
        # So keep all of the default values except plevs.
        self.granulate.remove("plevs")

    def check_values(self):
        plevs = self.plevs
        if not isinstance(plevs, list):
            msg = "plevs needs to be a list"
            raise RuntimeError(msg)

        if len(plevs) > 1:
            if monotonic(plevs):
                if plevs[0] > plevs[1]:
                    plevs = plevs[
                        ::-1
                    ]  # Force plevs to be monotonically increasing by reversing the list.
                    self.plevs = plevs
            else:
                msg = "plevs should be monotonically increasing or decreasing"
                raise RuntimeError(msg)
        else:
            msg = "At least 2 plevs needed"
            raise RuntimeError(msg)
