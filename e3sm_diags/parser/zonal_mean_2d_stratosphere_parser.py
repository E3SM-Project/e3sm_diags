from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)
from e3sm_diags.parser.core_parser import CoreParser


class ZonalMean2dStratosphereParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(  # type: ignore
            parameter_cls=ZonalMean2dStratosphereParameter, *args, **kwargs
        )
