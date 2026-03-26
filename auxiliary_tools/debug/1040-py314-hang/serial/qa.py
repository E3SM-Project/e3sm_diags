#!/usr/bin/env python3
"""Local Python driver for a complete model-vs-obs E3SM Diags run.

This replaces the original SLURM/bash wrapper with a direct Python entrypoint.
It stages year-filtered climatology files in a temporary local workspace and
then runs the same diagnostics via ``runner.run_diags()``.

Source: /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.serial.bash

Usage:
conda env create -f conda-dev/dev.yml -n ed_1040_py314
srun --pty --nodes=1 --time=02:00:00 /bin/bash
conda activate ed_1040_py314
python auxiliary_tools/debug/1040-py314-hang/serial/qa.py

Usage of source script with sbatch:
sbatch /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.serial.bash
"""
from __future__ import annotations

import inspect
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.mp_partition_parameter import MPpartitionParameter
from e3sm_diags.parameter.precip_pdf_parameter import PrecipPDFParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.tropical_subseasonal_parameter import (
    TropicalSubseasonalParameter,
)
from e3sm_diags.run import runner

os.environ["E3SM_DIAGS_DEBUG_HANG"] = "1"
os.environ["E3SM_DIAGS_ALBEDO_DEBUG"] = "1"
os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
os.environ["PYTHONFAULTHANDLER"] = "1"

CASE_NAME = "v3.LR.historical_0051"
SHORT_NAME = "v3.LR.historical_0051"
START_YEAR = 1985
END_YEAR = 2014
RUN_TYPE = "model_vs_obs"
RESULTS_SUBDIR = f"{RUN_TYPE}_{START_YEAR}-{END_YEAR}"

OBS_CLIMO_DIR = Path("/lcrc/group/e3sm/diagnostics/observations/Atm/climatology")
OBS_TS_DIR = Path("/lcrc/group/e3sm/diagnostics/observations/Atm/time-series")
OBS_TC_DIR = Path("/lcrc/group/e3sm/diagnostics/observations/Atm/tc-analysis")

ATM_ROOT = Path(
    "/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/"
    "v3.LR.historical_0051/post/atm/180x360_aave"
)
ROF_ROOT = Path(
    "/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/"
    "v3.LR.historical_0051/post/rof/native"
)
TC_ROOT = Path(
    "/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/"
    "v3.LR.historical_0051/post/atm"
)

SETS_TO_RUN = [
    "lat_lon",
]
_ALBEDO_DEBUG_MODE = os.environ.get("E3SM_DIAGS_ALBEDO_DEBUG", "1").lower() in (
    "1",
    "true",
    "yes",
    "on",
)

if not _ALBEDO_DEBUG_MODE:
    SETS_TO_RUN = [
        "lat_lon",
        "zonal_mean_xy",
        "zonal_mean_2d",
        "polar",
        "cosp_histogram",
        "meridional_mean_2d",
        "annual_cycle_zonal_mean",
        "enso_diags",
        "qbo",
        "diurnal_cycle",
        "zonal_mean_2d_stratosphere",
        "aerosol_aeronet",
        "mp_partition",
        "tropical_subseasonal",
        "precip_pdf",
        "tc_analysis",
        "streamflow",
    ]


SCRIPT_DIR = Path(__file__).resolve().parent
WORKDIR_ROOT: Path | None = None
RESULTS_DIR = SCRIPT_DIR / "output" / RESULTS_SUBDIR
NUM_WORKERS = 1
MULTIPROCESSING = False
KEEP_WORKDIR = False


def _maybe_wait_for_debugger():
    debugpy_port = os.environ.get("E3SM_DIAGS_DEBUGPY_PORT")
    if not debugpy_port:
        return

    import debugpy

    port = int(debugpy_port)
    debugpy.listen(("0.0.0.0", port))
    print(f"Waiting for debugger attach on 0.0.0.0:{port}")
    debugpy.wait_for_client()


def _print_runtime_paths():
    import e3sm_diags
    from e3sm_diags.driver.utils import dataset_xr

    print("Runtime module paths:")
    print(f"  e3sm_diags: {e3sm_diags.__file__}")
    print(f"  dataset_xr: {inspect.getsourcefile(dataset_xr.Dataset)}")


def main() -> int:
    start = time.time()
    results_dir = RESULTS_DIR.resolve()

    _maybe_wait_for_debugger()

    workdir = Path(
        tempfile.mkdtemp(
            prefix="tmp.e3sm_diags_atm_monthly_180x360_aave.",
            dir=WORKDIR_ROOT,
        )
    )
    try:
        run(
            workdir=workdir,
            results_dir=results_dir,
            multiprocessing=MULTIPROCESSING,
            num_workers=NUM_WORKERS,
        )
    finally:
        if KEEP_WORKDIR:
            print(f"Workspace retained at: {workdir}")
        else:
            shutil.rmtree(workdir, ignore_errors=True)

    elapsed = time.time() - start
    print(f"Completed local e3sm_diags run in {elapsed:.1f} seconds")
    print(f"Results written to: {results_dir}")
    return 0


def run(workdir: Path, results_dir: Path, multiprocessing: bool, num_workers: int) -> None:
    os.environ.setdefault("UVCDAT_ANONYMOUS_LOG", "False")
    results_dir.mkdir(parents=True, exist_ok=True)
    _print_runtime_paths()

    climo_dir, diurnal_dir = stage_inputs(workdir)
    params = build_params(
        climo_dir=climo_dir,
        diurnal_dir=diurnal_dir,
        results_dir=results_dir,
        multiprocessing=multiprocessing,
        num_workers=num_workers,
    )

    argv: list[str] = []
    argv.extend(["--no_viewer"])
    if multiprocessing:
        argv.extend(["--multiprocessing", "--num_workers", str(num_workers)])
    sys.argv.extend(argv)

    runner.sets_to_run = SETS_TO_RUN
    runner.run_diags(params)


def build_params(
    climo_dir: Path,
    diurnal_dir: Path,
    results_dir: Path,
    multiprocessing: bool,
    num_workers: int,
) -> list[object]:
    test_ts = ATM_ROOT / "ts" / "monthly" / "5yr"
    test_ts_daily = ATM_ROOT / "ts" / "daily" / "5yr"
    test_streamflow = ROF_ROOT / "ts" / "monthly" / "5yr"
    test_tc = TC_ROOT / f"tc-analysis_{START_YEAR}_{END_YEAR}"

    num_years = END_YEAR - START_YEAR + 1
    ref_start_yr = 1985

    param = CoreParameter()
    param.test_data_path = str(climo_dir)
    param.test_name = CASE_NAME
    param.short_test_name = SHORT_NAME
    param.reference_data_path = str(OBS_CLIMO_DIR)
    param.results_dir = str(results_dir)
    param.run_type = RUN_TYPE
    param.diff_title = "Model - Observations"
    param.output_format = ["png"]
    param.output_format_subplot = []
    param.multiprocessing = multiprocessing
    param.num_workers = num_workers
    param.no_viewer = True

    if _ALBEDO_DEBUG_MODE:
        param.sets = ["lat_lon"]
        param.variables = ["ALBEDO"]
        param.seasons = ["ANN", "DJF"]
        param.regions = ["global"]
        param.multiprocessing = False
        print("ALBEDO debug mode enabled:")
        print("  sets=['lat_lon']")
        print("  variables=['ALBEDO']")
        print("  seasons=['ANN', 'DJF']")
        print("  regions=['global']")
        print("  multiprocessing=False")
        return [param]

    params = [param]

    enso_param = EnsoDiagsParameter()
    enso_param.test_data_path = str(test_ts)
    enso_param.test_name = SHORT_NAME
    enso_param.test_start_yr = START_YEAR
    enso_param.test_end_yr = END_YEAR
    enso_param.reference_data_path = str(OBS_TS_DIR)
    enso_param.ref_start_yr = ref_start_yr
    enso_param.ref_end_yr = ref_start_yr + 10
    params.append(enso_param)

    trop_param = TropicalSubseasonalParameter()
    trop_param.test_data_path = str(test_ts_daily)
    trop_param.test_name = SHORT_NAME
    trop_param.test_start_yr = START_YEAR
    trop_param.test_end_yr = END_YEAR
    trop_param.reference_data_path = str(OBS_TS_DIR)
    trop_param.ref_start_yr = 2001
    trop_param.ref_end_yr = 2010
    params.append(trop_param)

    qbo_param = QboParameter()
    qbo_param.test_data_path = str(test_ts)
    qbo_param.test_name = SHORT_NAME
    qbo_param.test_start_yr = START_YEAR
    qbo_param.test_end_yr = END_YEAR
    qbo_param.ref_start_yr = ref_start_yr
    qbo_param.ref_end_yr = min(ref_start_yr + num_years - 1, 2014)
    qbo_param.reference_data_path = str(OBS_TS_DIR)
    params.append(qbo_param)

    mp_param = MPpartitionParameter()
    mp_param.test_data_path = str(test_ts)
    mp_param.test_name = SHORT_NAME
    mp_param.short_test_name = SHORT_NAME
    mp_param.test_start_yr = START_YEAR
    mp_param.test_end_yr = END_YEAR
    params.append(mp_param)

    precip_pdf_param = PrecipPDFParameter()
    precip_pdf_param.test_data_path = str(test_ts_daily)
    precip_pdf_param.test_name = SHORT_NAME
    precip_pdf_param.short_test_name = SHORT_NAME
    precip_pdf_param.test_start_yr = START_YEAR
    precip_pdf_param.test_end_yr = END_YEAR
    precip_pdf_param.reference_data_path = str(OBS_TS_DIR)
    params.append(precip_pdf_param)

    dc_param = DiurnalCycleParameter()
    dc_param.test_data_path = str(diurnal_dir)
    dc_param.short_test_name = SHORT_NAME
    dc_param.normalize_test_amp = False
    dc_param.reference_data_path = str(OBS_CLIMO_DIR)
    params.append(dc_param)

    streamflow_param = StreamflowParameter()
    streamflow_param.test_data_path = str(test_streamflow)
    streamflow_param.test_name = SHORT_NAME
    streamflow_param.test_start_yr = START_YEAR
    streamflow_param.test_end_yr = END_YEAR
    streamflow_param.reference_data_path = str(OBS_TS_DIR)
    streamflow_param.ref_start_yr = "1986"
    streamflow_param.ref_end_yr = "1995"
    params.append(streamflow_param)

    tc_param = TCAnalysisParameter()
    tc_param.test_data_path = str(test_tc)
    tc_param.short_test_name = SHORT_NAME
    tc_param.test_start_yr = str(START_YEAR)
    tc_param.test_end_yr = str(END_YEAR)
    tc_param.reference_data_path = str(OBS_TC_DIR)
    tc_param.ref_start_yr = "1979"
    tc_param.ref_end_yr = "2018"
    params.append(tc_param)

    return params


def stage_inputs(workdir: Path) -> tuple[Path, Path]:
    climo_dir = workdir / "climo"
    diurnal_dir = workdir / "climo_diurnal_8xdaily"

    year_tag = f"{START_YEAR}??_{END_YEAR}??"
    link_matches(
        ATM_ROOT / "clim" / "30yr",
        climo_dir,
        f"{CASE_NAME}_*_{year_tag}_climo.nc",
    )
    link_matches(
        ATM_ROOT / "clim_diurnal_8xdaily" / "30yr",
        diurnal_dir,
        f"{CASE_NAME}.*_*_{year_tag}_climo.nc",
    )

    return climo_dir, diurnal_dir


def link_matches(source_dir: Path, destination_dir: Path, pattern: str) -> None:
    matches = sorted(source_dir.glob(pattern))
    if not matches:
        raise FileNotFoundError(
            f"No files matched pattern {pattern!r} in {source_dir}"
        )

    destination_dir.mkdir(parents=True, exist_ok=True)
    for src in matches:
        target = destination_dir / src.name
        if target.exists() or target.is_symlink():
            target.unlink()
        target.symlink_to(src)


if __name__ == "__main__":
    raise SystemExit(main())
