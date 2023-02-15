from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter

from .core_parser_new import CoreParser


class TCAnalysisParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=TCAnalysisParameter, *args, **kwargs)  # type: ignore
