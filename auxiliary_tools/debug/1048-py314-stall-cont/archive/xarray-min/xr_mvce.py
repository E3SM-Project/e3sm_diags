#!/usr/bin/env python3
"""Minimal xarray reproduction case for Python 3.14 + xarray >=2026.01.0 stall.

Reproduces intermittent stall during multiprocessing open_mfdataset on NetCDF3 files.

Version boundary:
- xarray 2025.12.0: completes
- xarray >=2026.01.0: stalls in xarray/backends/locks.py
- xarray >=2026.01.0 with lock=False: completes

Expected behavior: All workers complete successfully
Observed behavior (xarray >=2026.01.0): Some runs stall, workers block in lock acquisition

Enhanced version: Each task opens 3 files (test climo + reference climo + mask file)
to match E3SM Diags workflow pattern that reproduces stall on NERSC/LCRC.

Usage:
  python xr_mvce.py                    # default lock (may stall on >=2026.01.0)
  USE_LOCK_FALSE=1 python xr_mvce.py   # lock=False workaround
"""
from __future__ import annotations

import multiprocessing as mp
import os
import sys
import time
from pathlib import Path
from typing import Any

import xarray as xr


TEST_CLIMO_FILE = "../v3.LR.historical_0051_ANN_198501_201412_climo.nc"

# On NERSC/LCRC: use actual obs climo file
# REF_CLIMO_FILE = "/global/cfs/projectdirs/e3sm/diagnostics/observations/Atm/climatology/CERES-EBAF-4.2_ANN_climo.nc"
# Fallback for local testing: reuse test file
REF_CLIMO_FILE = TEST_CLIMO_FILE

NUM_WORKERS = 8
NUM_ITERATIONS = 10
TIMEOUT_SECONDS = 120
USE_LOCK_FALSE = os.environ.get("USE_LOCK_FALSE", "0") == "1"

VARIABLES = ["ALBEDO", "TREFHT", "PRECT", "PSL", "QFLX", "SHFLX", "LHFLX", "TMQ"]


def worker_task(task_id: int) -> dict[str, Any]:
    """Simulate diagnostic task: open test + ref + mask files, access data, return."""
    pid = os.getpid()
    var_index = task_id % len(VARIABLES)
    var_name = VARIABLES[var_index]

    try:
        start = time.perf_counter()

        open_args = {
            "decode_times": True,
            "coords": "minimal",
            "compat": "override",
        }

        if USE_LOCK_FALSE:
            open_args["lock"] = False

        # Open test climo (mimics test_data_path open in E3SM Diags)
        ds_test = xr.open_mfdataset(paths=TEST_CLIMO_FILE, **open_args)

        # Open reference climo (mimics reference_data_path open in E3SM Diags)
        ds_ref = xr.open_mfdataset(paths=REF_CLIMO_FILE, **open_args)

        # Open mask file (mimics land/sea mask access in E3SM Diags)
        # Use open_dataset (not open_mfdataset) for single-file mask access
        mask_args = {"decode_times": True}
        if USE_LOCK_FALSE:
            mask_args["lock"] = False
        ds_mask = xr.open_dataset(TEST_CLIMO_FILE, **mask_args)

        # Load data from test file
        if var_name == "ALBEDO" and "FSNTOA" in ds_test and "SOLIN" in ds_test:
            _ = (ds_test["FSNTOA"] / ds_test["SOLIN"]).load()
        elif var_name in ds_test:
            _ = ds_test[var_name].load()
        else:
            first_var = list(ds_test.data_vars)[0]
            _ = ds_test[first_var].load()

        # Load data from reference file
        if var_name in ds_ref:
            _ = ds_ref[var_name].load()
        elif len(ds_ref.data_vars) > 0:
            first_var = list(ds_ref.data_vars)[0]
            _ = ds_ref[first_var].load()

        # Load mask data (LANDFRAC or fallback)
        if "LANDFRAC" in ds_mask:
            _ = ds_mask["LANDFRAC"].load()
        elif "OCNFRAC" in ds_mask:
            _ = ds_mask["OCNFRAC"].load()
        elif len(ds_mask.data_vars) > 0:
            first_var = list(ds_mask.data_vars)[0]
            _ = ds_mask[first_var].load()

        ds_test.close()
        ds_ref.close()
        ds_mask.close()

        elapsed = time.perf_counter() - start

        return {
            "status": "ok",
            "task_id": task_id,
            "pid": pid,
            "var": var_name,
            "elapsed": round(elapsed, 3),
        }

    except Exception as exc:
        return {
            "status": "error",
            "task_id": task_id,
            "pid": pid,
            "var": var_name,
            "error": repr(exc),
        }


def main() -> int:
    test_path = Path(TEST_CLIMO_FILE)
    ref_path = Path(REF_CLIMO_FILE)

    if not test_path.exists():
        print(f"Error: Test climo file not found: {test_path.resolve()}")
        return 1

    if not ref_path.exists():
        print(f"Error: Ref climo file not found: {ref_path.resolve()}")
        return 1

    print(f"xarray version: {xr.__version__}")
    print(f"Python version: {sys.version}")
    print(f"Test climo file: {test_path.resolve()}")
    print(f"Ref climo file: {ref_path.resolve()}")
    print(f"Workers: {NUM_WORKERS}")
    print(f"Tasks per worker: {NUM_ITERATIONS}")
    print(f"Total tasks: {NUM_WORKERS * NUM_ITERATIONS}")
    print(f"Files opened per task: 3 (test + ref + mask)")
    print(f"lock=False: {USE_LOCK_FALSE}")
    print(f"Timeout: {TIMEOUT_SECONDS}s")
    print()

    task_ids = list(range(NUM_WORKERS * NUM_ITERATIONS))

    ctx = mp.get_context("fork")

    start = time.perf_counter()
    with ctx.Pool(processes=NUM_WORKERS) as pool:
        try:
            results = pool.map(worker_task, task_ids, chunksize=1)
        except mp.TimeoutError:
            print("\n❌ TIMEOUT: Workers did not complete within timeout period")
            return 1

    elapsed = time.perf_counter() - start

    ok_results = [r for r in results if r.get("status") == "ok"]
    error_results = [r for r in results if r.get("status") == "error"]

    print(f"✅ Completed: {len(ok_results)}/{len(task_ids)} tasks")
    print(f"❌ Errors: {len(error_results)}")
    print(f"⏱️  Total time: {elapsed:.1f}s")

    if error_results:
        print("\nErrors:")
        for err in error_results[:5]:
            print(f"  Task {err['task_id']} ({err['var']}): {err['error']}")

    if len(ok_results) < len(task_ids):
        print(f"\n❌ INCOMPLETE: Only {len(ok_results)}/{len(task_ids)} tasks completed")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
