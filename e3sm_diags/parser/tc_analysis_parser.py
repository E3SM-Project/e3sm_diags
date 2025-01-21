from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parser.core_parser import CoreParser


class TCAnalysisParser(CoreParser):
    def __init__(self, *args, **kwargs):
        # FIXME: B026 Star-arg unpacking after a keyword argument is strongly discouraged
        super().__init__(parameter_cls=TCAnalysisParameter, *args, **kwargs)  # type: ignore # noqa: B026
