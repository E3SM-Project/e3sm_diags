#!/usr/bin/env python3
"""Local Python driver for a complete model-vs-obs E3SM Diags run.

This replaces the original SLURM/bash wrapper with a direct Python entrypoint.
It stages year-filtered climatology files in a temporary local workspace and
then runs the same diagnostics via ``runner.run_diags()``.

Source: /lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.serial.bash

Usage:
conda env create -f conda-dev/dev.yml -n <ed_main_20e72b2e | ed_main_f80253a0>
srun --pty --nodes=1 --time=02:00:00 /bin/bash
conda activate <ed_main_20e72b2e | ed_main_f80253a0>
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
from e3sm_diags.run import runner


CASE_NAME = "v3.LR.historical_0051"
SHORT_NAME = "v3.LR.historical_0051"
START_YEAR = 1985
END_YEAR = 2014
RUN_TYPE = "model_vs_obs"
RESULTS_SUBDIR = f"{RUN_TYPE}_{START_YEAR}-{END_YEAR}"

OBS_CLIMO_DIR = Path("/global/cfs/projectdirs/e3sm/diagnosticsAtm/climatology")
OBS_TS_DIR = Path("/global/cfs/projectdirs/e3sm/diagnosticsAtm/time-series")
OBS_TC_DIR = Path("/global/cfs/projectdirs/e3sm/diagnosticsAtm/tc-analysis")

ATM_ROOT = Path("/global/cfs/projectdirs/e3sm/vo13/1048-py314-stall-cont/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave")

SETS_TO_RUN = ["lat_lon"]
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
CONDA_ENV = os.environ.get("CONDA_DEFAULT_ENV", "unknown_env")
RESULTS_DIR = SCRIPT_DIR / "output" / CONDA_ENV / RESULTS_SUBDIR
NUM_WORKERS = 4
MULTIPROCESSING = True
KEEP_WORKDIR = False


def _print_runtime_paths():
    import e3sm_diags
    from e3sm_diags.driver.utils import dataset_xr

    print("Runtime module paths:")
    print(f"  e3sm_diags: {e3sm_diags.__file__}")
    print(f"  dataset_xr: {inspect.getsourcefile(dataset_xr.Dataset)}")


def main() -> int:
    start = time.time()
    results_dir = RESULTS_DIR.resolve()


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

    climo_dir = stage_inputs(workdir)
    params = build_params(
        climo_dir=climo_dir,
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
    results_dir: Path,
    multiprocessing: bool,
    num_workers: int,
) -> list[object]:

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

    param.sets = ["lat_lon"]
    param.variables = ["ALBEDO"]
    param.seasons = ["ANN", "DJF"]
    param.regions = ["global"]
    param.multiprocessing = True
    print("ALBEDO debug mode enabled:")
    print("  sets=['lat_lon']")
    print("  variables=['ALBEDO']")
    print("  seasons=['ANN', 'DJF']")
    print("  regions=['global']")
    print("  multiprocessing=True")
    return [param]



def stage_inputs(workdir: Path) -> tuple[Path, Path]:
    climo_dir = workdir / "climo"

    year_tag = f"{START_YEAR}??_{END_YEAR}??"
    link_matches(
        ATM_ROOT / "clim" / "30yr",
        climo_dir,
        f"{CASE_NAME}_*_{year_tag}_climo.nc",
    )

    return climo_dir


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
