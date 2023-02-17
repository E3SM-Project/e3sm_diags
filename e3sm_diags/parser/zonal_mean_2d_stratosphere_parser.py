from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)
from e3sm_diags.parser.zonal_mean_2d_parser import ZonalMean2dParser


class ZonalMean2dStratosphereParser(ZonalMean2dParser):
    def __init__(self, *args, **kwargs):
        # `parameter_cls` must be last because ZonalMean2dParser already
        # defines `parameter_cls=ZonalMean2dParameter`, so we need to override
        # it.
        super().__init__(
            *args,
            **kwargs,
            parameter_cls=ZonalMean2dStratosphereParameter,
        )
