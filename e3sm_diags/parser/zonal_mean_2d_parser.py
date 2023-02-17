from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.parser.core_parser import CoreParser


class ZonalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=ZonalMean2dParameter, *args, **kwargs)  # type: ignore
