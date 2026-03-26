#!/usr/bin/env python3
"""Fast polar-only forkserver repro for the Python 3.14 hang investigation.

This replaces the longer model-vs-obs wrapper with a focused local entrypoint
that stages year-filtered climatology files in a temporary workspace and runs a
small `polar` diagnostic repro via ``runner.run_diags()``.

Default repro:
- set: ``polar``
- variable: ``T``
- seasons: ``ANN,DJF``
- regions: ``polar_S,polar_N``
- multiprocessing: enabled
- num_workers: ``4``

Useful overrides:
- ``E3SM_DIAGS_REPRO_VAR=Z3``
- ``E3SM_DIAGS_REPRO_SEASONS=ANN,DJF``
- ``E3SM_DIAGS_REPRO_REGIONS=polar_S,polar_N``
- ``E3SM_DIAGS_REPRO_PLEVS=500.0``
- ``E3SM_DIAGS_REPRO_WORKERS=8``

Usage:
conda env create -f conda-dev/dev.yml -n ed_1040_py314
srun --pty --nodes=1 --time=02:00:00 /bin/bash
conda activate ed_1040_py314
python auxiliary_tools/debug/1040-py314-hang/parallel/qa.py

Reference script:
/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.py314_tom_branch.bash
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

SETS_TO_RUN = ["polar"]
REPRO_VAR = "T"
REPRO_SEASONS = ["ANN"]
REPRO_REGIONS = ["polar_S","polar_N"]
REPRO_PLEVS_RAW = "200.0"
SCRIPT_DIR = Path(__file__).resolve().parent
WORKDIR_ROOT: Path | None = None
RESULTS_DIR = SCRIPT_DIR / "output" / RESULTS_SUBDIR
KEEP_WORKDIR = False

# NOTE: Update these variables for testing different cases.
MULTIPROCESSING = True
NUM_WORKERS = 8



def _parse_plevs(raw: str) -> list[float]:
    if raw == "":
        return []

    return [float(item.strip()) for item in raw.split(",") if item.strip()]


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
    print("Quick repro configuration:")
    print(f"  sets={SETS_TO_RUN}")
    print(f"  variable={REPRO_VAR}")
    print(f"  seasons={REPRO_SEASONS}")
    print(f"  regions={REPRO_REGIONS}")
    print(f"  plevs={_parse_plevs(REPRO_PLEVS_RAW)}")
    print(f"  multiprocessing={MULTIPROCESSING}")
    print(f"  num_workers={NUM_WORKERS}")


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
    param.sets = ["polar"]
    param.variables = [REPRO_VAR]
    param.seasons = REPRO_SEASONS
    param.regions = REPRO_REGIONS
    param.plevs = _parse_plevs(REPRO_PLEVS_RAW)
    param.no_viewer = True

    return [param]


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
