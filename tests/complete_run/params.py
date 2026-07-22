"""Build parameter objects for the manual complete-run workflow.

This module centralizes the filesystem defaults and parameter-construction
logic used by the manual complete-run diagnostics flow. It exists so the
current branch's working parameter setup can be preserved in one place and
reused by the run entrypoint without hardcoding behavior into top-level
scripts.
"""

from __future__ import annotations

from dataclasses import dataclass

from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.tropical_subseasonal_parameter import (
    TropicalSubseasonalParameter,
)
from tests.complete_run.helpers import append_run_suffix


@dataclass(frozen=True)
class CompleteRunPaths:
    """Filesystem inputs for a manual complete-run workflow."""

    results_dir: str
    test_climo: str
    test_ts: str
    test_ts_daily_dir: str
    test_diurnal_climo: str
    test_streamflow_ts: str
    test_tc_analysis: str
    test_arm_site: str
    ref_climo: str
    ref_ts: str
    ref_tc_analysis: str
    ref_arm: str


@dataclass(frozen=True)
class CompleteRunConfig:
    """Runtime configuration for a manual complete run."""

    case: str
    short_name: str
    start_yr: str
    end_yr: str
    num_workers: int
    save_netcdf: bool
    paths: CompleteRunPaths


DEFAULT_CASE = "extendedOutput.v3.LR.historical_0101"
DEFAULT_SHORT_NAME = "v3.LR.historical_0101"
DEFAULT_START_YEAR = "2000"
DEFAULT_END_YEAR = "2014"
DEFAULT_NUM_WORKERS = 24

# Default Paths for input data and results
# ----------------------------------------
# TODO: Do we still want to keep using this directory for test input data?
DEFAULT_TEST_INPUT_PATH = (
    f"/global/cfs/cdirs/e3sm/chengzhu/tutorial2024/{DEFAULT_SHORT_NAME}"
)
DEFAULT_REF_INPUT_PATH = "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm"
DEFAULT_RESULTS_DIR = "/global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test"


def build_default_paths() -> CompleteRunPaths:
    return CompleteRunPaths(
        results_dir=append_run_suffix(DEFAULT_RESULTS_DIR),
        test_climo=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/180x360_aave/clim/15yr",
        test_ts=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/180x360_aave/ts/monthly/15yr",
        test_ts_daily_dir=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/180x360_aave/ts/daily/15yr",
        test_diurnal_climo=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/180x360_aave/clim_diurnal_8xdaily/",
        test_streamflow_ts=f"{DEFAULT_TEST_INPUT_PATH}/post/rof/native/ts/monthly/15yr/",
        test_tc_analysis=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/tc-analysis_2000_2014",
        test_arm_site=f"{DEFAULT_TEST_INPUT_PATH}/post/atm/site",
        ref_climo=f"{DEFAULT_REF_INPUT_PATH}/climatology/",
        ref_ts=f"{DEFAULT_REF_INPUT_PATH}/time-series",
        ref_tc_analysis=f"{DEFAULT_REF_INPUT_PATH}/tc-analysis/",
        ref_arm=f"{DEFAULT_REF_INPUT_PATH}/arm-diags-data",
    )


DEFAULT_PATHS = build_default_paths()


def build_complete_run_config(
    *,
    case: str = DEFAULT_CASE,
    short_name: str = DEFAULT_SHORT_NAME,
    start_yr: str = DEFAULT_START_YEAR,
    end_yr: str = DEFAULT_END_YEAR,
    num_workers: int = DEFAULT_NUM_WORKERS,
    save_netcdf: bool = True,
    paths: CompleteRunPaths | None = None,
) -> CompleteRunConfig:
    """Build a validated runtime configuration for a complete run.

    Parameters
    ----------
    case : str, optional
        Test-case name passed to the core parameter.
    short_name : str, optional
        Short display name for the test case.
    start_yr : str, optional
        Test data start year for time-series diagnostics.
    end_yr : str, optional
        Test data end year for time-series diagnostics.
    num_workers : int, optional
        Worker count for multiprocessing diagnostics.
    save_netcdf : bool, optional
        Whether to write netCDF outputs for the manual comparison workflow.
    paths : CompleteRunPaths, optional
        Input and output paths used by the diagnostics.

    Returns
    -------
    CompleteRunConfig
        The runtime configuration.

    Raises
    ------
    ValueError
        If ``num_workers`` is less than 1.
    """
    if num_workers < 1:
        raise ValueError("num_workers must be greater than or equal to 1.")

    if paths is None:
        paths = build_default_paths()

    return CompleteRunConfig(
        case=case,
        short_name=short_name,
        start_yr=start_yr,
        end_yr=end_yr,
        num_workers=num_workers,
        save_netcdf=save_netcdf,
        paths=paths,
    )


def build_complete_run_params(config: CompleteRunConfig) -> list[CoreParameter]:
    """Build parameter objects for the manual complete-run workflow.

    Parameters
    ----------
    config : CompleteRunConfig
        Runtime settings and filesystem locations for the complete run.

    Returns
    -------
    list[CoreParameter]
        The parameter objects used by ``runner.run_diags()``.
    """
    paths = config.paths
    params: list[CoreParameter] = []

    param = CoreParameter()
    param.test_data_path = paths.test_climo
    param.test_name = config.case
    param.short_test_name = config.short_name
    param.reference_data_path = paths.ref_climo
    param.results_dir = paths.results_dir
    param.run_type = "model_vs_obs"
    param.diff_title = "Model - Observations"
    param.output_format = ["png"]
    param.output_format_subplot = []
    param.multiprocessing = config.num_workers > 1
    param.num_workers = config.num_workers
    param.save_netcdf = config.save_netcdf
    param.seasons = ["ANN"]
    params.append(param)

    enso_param = EnsoDiagsParameter()
    enso_param.test_data_path = paths.test_ts
    enso_param.test_start_yr = config.start_yr
    enso_param.test_end_yr = config.end_yr
    enso_param.reference_data_path = paths.ref_ts
    enso_param.ref_start_yr = config.start_yr
    enso_param.ref_end_yr = config.end_yr
    enso_param.save_netcdf = config.save_netcdf
    params.append(enso_param)

    trop_param = TropicalSubseasonalParameter()
    trop_param.test_data_path = paths.test_ts_daily_dir
    trop_param.test_start_yr = config.start_yr
    trop_param.test_end_yr = config.end_yr
    trop_param.reference_data_path = paths.ref_ts
    trop_param.ref_start_yr = "2001"
    trop_param.ref_end_yr = "2010"
    trop_param.save_netcdf = config.save_netcdf
    params.append(trop_param)

    qbo_param = QboParameter()
    qbo_param.test_data_path = paths.test_ts
    qbo_param.test_start_yr = config.start_yr
    qbo_param.test_end_yr = config.end_yr
    qbo_param.ref_start_yr = config.start_yr
    qbo_param.ref_end_yr = config.end_yr
    qbo_param.reference_data_path = paths.ref_ts
    qbo_param.save_netcdf = config.save_netcdf
    params.append(qbo_param)

    dc_param = DiurnalCycleParameter()
    dc_param.test_data_path = paths.test_diurnal_climo
    dc_param.normalize_test_amp = False
    dc_param.reference_data_path = paths.ref_climo
    dc_param.save_netcdf = config.save_netcdf
    params.append(dc_param)

    streamflow_param = StreamflowParameter()
    streamflow_param.test_data_path = paths.test_streamflow_ts
    streamflow_param.test_start_yr = config.start_yr
    streamflow_param.test_end_yr = config.end_yr
    streamflow_param.reference_data_path = paths.ref_ts
    streamflow_param.ref_start_yr = "1986"
    streamflow_param.ref_end_yr = "1995"
    streamflow_param.save_netcdf = config.save_netcdf
    params.append(streamflow_param)

    tc_param = TCAnalysisParameter()
    tc_param.test_data_path = paths.test_tc_analysis
    tc_param.test_start_yr = config.start_yr
    tc_param.test_end_yr = config.end_yr
    tc_param.reference_data_path = paths.ref_tc_analysis
    tc_param.ref_start_yr = "1979"
    tc_param.ref_end_yr = "2018"
    tc_param.save_netcdf = config.save_netcdf
    params.append(tc_param)

    arm_param = ARMDiagsParameter()
    arm_param.reference_data_path = paths.ref_arm
    arm_param.ref_name = "armdiags"
    arm_param.test_data_path = paths.test_arm_site
    arm_param.test_name = config.short_name
    arm_param.test_start_yr = config.start_yr
    arm_param.test_end_yr = config.end_yr
    arm_param.ref_start_yr = "0001"
    arm_param.ref_end_yr = "0001"
    arm_param.save_netcdf = config.save_netcdf
    params.append(arm_param)

    return params
