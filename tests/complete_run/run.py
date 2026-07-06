"""Run the manual complete E3SM diagnostics workflow.

This module provides the manual command-line entrypoint for launching the
complete-run diagnostics flow against HPC-resident input data. It wraps the
current complete-run parameter and set configuration in a safe interface so
developers can reproduce broad diagnostics coverage without import-time side
effects, ad hoc script edits, or pytest collection.

Usage
-----
On NERSC:

    # Request an interactive node with a 4-hour walltime. Adjust the account,
    salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account=e3sm

    # Activate the conda environment used for E3SM Diags development.
    conda activate <e3sm_diags_env>

    # Run the complete-run workflow with default parameters and sets.
    # Use --help for more details on the available flags and their usage.
    python -m tests.complete_run.run
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Sequence

from e3sm_diags.logger import _setup_child_logger, _setup_root_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner
from tests.complete_run.params import (
    DEFAULT_CASE,
    DEFAULT_END_YEAR,
    DEFAULT_NUM_WORKERS,
    DEFAULT_SHORT_NAME,
    DEFAULT_START_YEAR,
    CompleteRunPaths,
    build_complete_run_config,
    build_complete_run_params,
    build_default_paths,
)

logger = _setup_child_logger(__name__)

DEFAULT_SETS_TO_RUN = [
    "lat_lon",
    "zonal_mean_xy",
    "zonal_mean_2d",
    "polar",
    "cosp_histogram",
    "meridional_mean_2d",
    "enso_diags",
    "qbo",
    "diurnal_cycle",
    "annual_cycle_zonal_mean",
    "streamflow",
    "zonal_mean_2d_stratosphere",
    "arm_diags",
    "tc_analysis",
    "aerosol_aeronet",
    "aerosol_budget",
    "tropical_subseasonal",
]


def main(argv: Sequence[str] | None = None) -> int:
    """CLI entrypoint for the manual complete-run workflow."""
    _setup_root_logger()

    parser = _build_parser()
    args = parser.parse_args(argv)
    results = _run_complete_run(args)

    if results:
        logger.info("Complete-run finished with %s parameter results.", len(results))
    else:
        logger.warning(
            "Complete-run finished without explicit parameter results. "
            "Check the run log under the results directory if needed."
        )

    return 0


def _build_parser() -> argparse.ArgumentParser:
    """Build the CLI parser for the manual complete-run workflow."""
    default_paths = build_default_paths()

    parser = argparse.ArgumentParser(
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Run the manual complete E3SM diagnostics workflow used for HPC "
            "validation and netCDF result generation."
        ),
        epilog=(
            "Example: python -m tests.complete_run.run "
            f"--results-dir {default_paths.results_dir} "
            "--set lat_lon --set enso_diags"
        ),
    )
    parser.add_argument("--case", default=DEFAULT_CASE, help="Test-case name.")
    parser.add_argument(
        "--short-name",
        default=DEFAULT_SHORT_NAME,
        help="Short display name for the test case.",
    )
    parser.add_argument(
        "--results-dir",
        default=default_paths.results_dir,
        help="Root output directory for this complete run.",
    )
    parser.add_argument(
        "--test-climo",
        default=default_paths.test_climo,
        help="Test climatology input directory.",
    )
    parser.add_argument(
        "--test-ts",
        default=default_paths.test_ts,
        help="Test monthly time-series input directory.",
    )
    parser.add_argument(
        "--test-ts-daily-dir",
        default=default_paths.test_ts_daily_dir,
        help="Test daily time-series input directory.",
    )
    parser.add_argument(
        "--test-diurnal-climo",
        default=default_paths.test_diurnal_climo,
        help="Test diurnal climatology input directory.",
    )
    parser.add_argument(
        "--test-streamflow-ts",
        default=default_paths.test_streamflow_ts,
        help="Test streamflow time-series input directory.",
    )
    parser.add_argument(
        "--test-tc-analysis",
        default=default_paths.test_tc_analysis,
        help="Test tropical cyclone analysis input directory.",
    )
    parser.add_argument(
        "--test-arm-site",
        default=default_paths.test_arm_site,
        help="Test ARM site input directory.",
    )
    parser.add_argument(
        "--ref-climo",
        default=default_paths.ref_climo,
        help="Reference climatology input directory.",
    )
    parser.add_argument(
        "--ref-ts",
        default=default_paths.ref_ts,
        help="Reference time-series input directory.",
    )
    parser.add_argument(
        "--ref-tc-analysis",
        default=default_paths.ref_tc_analysis,
        help="Reference tropical cyclone analysis input directory.",
    )
    parser.add_argument(
        "--ref-arm",
        default=default_paths.ref_arm,
        help="Reference ARM diagnostics input directory.",
    )
    parser.add_argument(
        "--start-yr",
        default=DEFAULT_START_YEAR,
        help="Start year for test time-series diagnostics.",
    )
    parser.add_argument(
        "--end-yr",
        default=DEFAULT_END_YEAR,
        help="End year for test time-series diagnostics.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=DEFAULT_NUM_WORKERS,
        help="Worker count for multiprocessing diagnostics.",
    )
    parser.add_argument(
        "--set",
        dest="sets_to_run",
        action="append",
        choices=DEFAULT_SETS_TO_RUN,
        help="Optional subset of diagnostic sets to run. Repeat to select multiple.",
    )
    parser.add_argument(
        "--save-netcdf",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Write netCDF outputs for later comparison.",
    )

    return parser


def _build_paths_from_args(args: argparse.Namespace) -> CompleteRunPaths:
    return CompleteRunPaths(
        results_dir=args.results_dir,
        test_climo=args.test_climo,
        test_ts=args.test_ts,
        test_ts_daily_dir=args.test_ts_daily_dir,
        test_diurnal_climo=args.test_diurnal_climo,
        test_streamflow_ts=args.test_streamflow_ts,
        test_tc_analysis=args.test_tc_analysis,
        test_arm_site=args.test_arm_site,
        ref_climo=args.ref_climo,
        ref_ts=args.ref_ts,
        ref_tc_analysis=args.ref_tc_analysis,
        ref_arm=args.ref_arm,
    )


def _validate_input_paths(paths: CompleteRunPaths) -> None:
    """Validate that required manual input directories exist.

    Parameters
    ----------
    paths : CompleteRunPaths
        Manual input and output path configuration.

    Raises
    ------
    FileNotFoundError
        If any required input directory does not exist.
    """
    path_map = {
        "test_climo": paths.test_climo,
        "test_ts": paths.test_ts,
        "test_ts_daily_dir": paths.test_ts_daily_dir,
        "test_diurnal_climo": paths.test_diurnal_climo,
        "test_streamflow_ts": paths.test_streamflow_ts,
        "test_tc_analysis": paths.test_tc_analysis,
        "test_arm_site": paths.test_arm_site,
        "ref_climo": paths.ref_climo,
        "ref_ts": paths.ref_ts,
        "ref_tc_analysis": paths.ref_tc_analysis,
        "ref_arm": paths.ref_arm,
    }

    missing_inputs = [
        f"{name}={path_str}"
        for name, path_str in path_map.items()
        if not Path(path_str).exists()
    ]
    if missing_inputs:
        missing_message = "\n".join(missing_inputs)
        raise FileNotFoundError(
            f"One or more complete-run input paths do not exist:\n{missing_message}"
        )


def _run_complete_run(args: argparse.Namespace) -> list[CoreParameter] | None:
    """Run the configured diagnostics workflow.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed CLI arguments.

    Returns
    -------
    list | None
        The result returned by ``runner.run_diags()``.
    """
    config = build_complete_run_config(
        case=args.case,
        short_name=args.short_name,
        start_yr=args.start_yr,
        end_yr=args.end_yr,
        num_workers=args.num_workers,
        save_netcdf=args.save_netcdf,
        paths=_build_paths_from_args(args),
    )
    _validate_input_paths(config.paths)
    params = build_complete_run_params(config)
    selected_sets = args.sets_to_run or DEFAULT_SETS_TO_RUN

    logger.info("Running complete-run workflow.")
    logger.info("Results dir: %s", config.paths.results_dir)
    logger.info("Selected sets (%s): %s", len(selected_sets), ", ".join(selected_sets))
    logger.info(
        "Years: test %s-%s | workers: %s | save_netcdf: %s",
        config.start_yr,
        config.end_yr,
        config.num_workers,
        config.save_netcdf,
    )

    runner.sets_to_run = selected_sets
    return runner.run_diags(params)


if __name__ == "__main__":
    raise SystemExit(main())
