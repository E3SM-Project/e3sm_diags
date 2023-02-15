from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter

from .core_parser_new import CoreParser


class ZonalMean2dParser(CoreParser):
    def __init__(self, *args, **kwargs):
        super().__init__(parameter_cls=ZonalMean2dParameter, *args, **kwargs)  # type: ignore

    def load_default_args(self):
        # The parameters unique to ZonalMean2dParameter are added here.
        self.add_argument(
            "--plevs",
            type=float,
            nargs="+",
            dest="plevs",
            help="Selected pressure level.[take list as input]",
            required=False,
        )

        self.add_argument(
            "--plot_plevs",
            dest="plot_plevs",
            help="plot specified plevs",
            action="store_const",
            const=True,
            required=False,
        )

        self.add_argument(
            "--plot_log_plevs",
            dest="plot_log_plevs",
            help="plot plevs on log-scale",
            action="store_const",
            const=True,
            required=False,
        )
