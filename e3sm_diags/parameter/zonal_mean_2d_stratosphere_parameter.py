import copy

import numpy as np

from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter

DEFAULT_PLEVS = np.logspace(0, 2.0, num=10).tolist()


class ZonalMean2dStratosphereParameter(ZonalMean2dParameter):
    def __init__(self):
        super(ZonalMean2dStratosphereParameter, self).__init__()
        # Override existing attributes
        # =============================
        self.plevs = copy.deepcopy(DEFAULT_PLEVS)  # type: ignore
        self.plot_log_plevs = True
