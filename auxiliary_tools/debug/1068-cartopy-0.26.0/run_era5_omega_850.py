"""Render only the complete-run ERA5 OMEGA 850 hPa ANN lat-lon diagnostic.

This reproduces the inputs and plot settings for
``ERA5-OMEGA-850-ANN-global.png`` from the Cartopy 0.26 complete-run comparison.
Use the same output directory with two environments to isolate rendering changes:

    python run_era5_omega_850.py --results-dir results-cartopy-026

The script reads the shared complete-run inputs but writes only to the selected
results directory. Run it from an allocated NERSC compute node.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from e3sm_diags.logger import _setup_root_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner


TEST_DATA_PATH = (
    "/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/v3.LR.historical_0101/"
    "post/atm/180x360_aave/clim/15yr"
)
REFERENCE_DATA_PATH = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology"
DEFAULT_RESULTS_DIR = Path(__file__).with_name("results-era5-omega-850")

CONTOUR_LEVELS = [
    -220.0,
    -180.0,
    -140.0,
    -100.0,
    -60.0,
    -20.0,
    -10.0,
    10.0,
    20.0,
    60.0,
    100.0,
    140.0,
    180.0,
    220.0,
]
DIFF_LEVELS = [
    -80.0,
    -60.0,
    -40.0,
    -32.0,
    -24.0,
    -16.0,
    -8.0,
    -4.0,
    4.0,
    8.0,
    16.0,
    24.0,
    32.0,
    40.0,
    60.0,
    80.0,
]
DIFF_LEVELS_ZERO_BOUNDARY = [
    -80.0,
    -60.0,
    -40.0,
    -32.0,
    -24.0,
    -16.0,
    -8.0,
    -4.0,
    0.0,
    4.0,
    8.0,
    16.0,
    24.0,
    32.0,
    40.0,
    60.0,
    80.0,
]
DIFF_LEVELS_COARSE = [
    -80.0,
    -60.0,
    -40.0,
    -20.0,
    -10.0,
    -5.0,
    0.0,
    5.0,
    10.0,
    20.0,
    40.0,
    60.0,
    80.0,
]
DIFF_LEVEL_EXPERIMENTS = {
    "standard": DIFF_LEVELS,
    "zero-boundary": DIFF_LEVELS_ZERO_BOUNDARY,
    "coarse": DIFF_LEVELS_COARSE,
}


def _build_parser() -> argparse.ArgumentParser:
    """Build the command-line parser."""
    parser = argparse.ArgumentParser(description=__doc__)
    # Jupyter starts kernels with ``-f``/``--f=<connection-file>``. Accept the
    # launcher argument so this script can be run from an interactive console.
    parser.add_argument(
        "-f", "--f", dest="_jupyter_connection_file", help=argparse.SUPPRESS
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=DEFAULT_RESULTS_DIR,
        help="Directory in which the single diagnostic result is written.",
    )
    parser.add_argument(
        "--diff-level-experiment",
        choices=[*DIFF_LEVEL_EXPERIMENTS, "all"],
        default="standard",
        help=(
            "Difference contour-level set to render. Use 'all' to render every "
            "experiment in named subdirectories under --results-dir."
        ),
    )
    return parser


def build_parameter(results_dir: Path, diff_levels: list[float]) -> CoreParameter:
    """Build the parameter matching the complete-run ERA5 OMEGA 850 hPa plot."""
    parameter = CoreParameter()
    parameter.test_data_path = TEST_DATA_PATH
    parameter.reference_data_path = REFERENCE_DATA_PATH
    parameter.results_dir = str(results_dir)
    parameter.run_type = "model_vs_obs"
    parameter.case_id = "ERA5"
    parameter.test_name = "extendedOutput.v3.LR.historical_0101"
    parameter.short_test_name = "v3.LR.historical_0101"
    parameter.ref_name = "ERA5"
    parameter.reference_name = "ERA5 Reanalysis"
    parameter.variables = ["OMEGA"]
    parameter.seasons = ["ANN"]
    parameter.plevs = [850.0]
    parameter.regions = ["global"]
    parameter.test_colormap = "PiYG_r"
    parameter.reference_colormap = "PiYG_r"
    # Set explicitly rather than relying on default-config resolution. The
    # complete-run image uses this blue-white-red difference palette.
    parameter.diff_colormap = "diverging_bwr.rgb"
    parameter.contour_levels = CONTOUR_LEVELS
    parameter.diff_levels = diff_levels
    parameter.regrid_method = "bilinear"
    parameter.diff_title = "Model - Observations"
    parameter.output_format = ["png"]
    parameter.output_format_subplot = []
    parameter.save_netcdf = False
    parameter.no_viewer = True
    parameter.multiprocessing = False
    return parameter


def main(argv: Sequence[str] | None = None) -> int:
    """Run the isolated lat-lon diagnostic."""
    _setup_root_logger()
    args = _build_parser().parse_args(argv)
    runner.sets_to_run = ["lat_lon"]

    experiments = DIFF_LEVEL_EXPERIMENTS.items()
    if args.diff_level_experiment != "all":
        experiments = [
            (
                args.diff_level_experiment,
                DIFF_LEVEL_EXPERIMENTS[args.diff_level_experiment],
            )
        ]

    for experiment_name, diff_levels in experiments:
        results_dir = args.results_dir
        if args.diff_level_experiment == "all":
            results_dir /= experiment_name
        runner.run_diags([build_parameter(results_dir, diff_levels)])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
