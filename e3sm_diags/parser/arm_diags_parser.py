from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter

from .core_parser_new import CoreParser


class ARMDiagsParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=ARMDiagsParameter, *args, **kwargs)  # type: ignore

    def load_default_args(self):
        self.add_argument(
            "--ref_names",
            type=str,
            nargs="+",
            dest="ref_names",
            help="List of reference names.",
            required=False,
        )
