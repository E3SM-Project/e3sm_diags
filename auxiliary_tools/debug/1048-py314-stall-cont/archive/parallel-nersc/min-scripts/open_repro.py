#!/usr/bin/env python3
from __future__ import annotations

import json
import multiprocessing as mp
import os
import queue
import time
from pathlib import Path
from typing import Any

import xarray as xr
import xcdat as xc


ANN_CLIMO_FILE = Path(
    "/tmp/tmp.e3sm_diags_atm_monthly_180x360_aave.lcui1k_l/climo/"
    "v3.LR.historical_0051_ANN_198501_201412_climo.nc"
)
NUM_WORKERS = int(os.environ.get("OPEN_REPRO_WORKERS", "8"))
NUM_ITERATIONS = int(os.environ.get("OPEN_REPRO_ITERATIONS", "20"))
WORKER_TIMEOUT_SECONDS = int(os.environ.get("OPEN_REPRO_TIMEOUT", "45"))


def _open_with_xc_mfdataset(path: Path) -> dict[str, Any]:
    ds = xc.open_mfdataset(
        paths=str(path),
        decode_times=True,
        add_bounds=["X", "Y", "Z"],
        coords="minimal",
        compat="override",
    )
    try:
        return {
            "dims": dict(ds.sizes),
            "data_vars": sorted(str(key) for key in ds.data_vars.keys())[:5],
        }
    finally:
        ds.close()


def _open_with_xc_dataset(path: Path) -> dict[str, Any]:
    ds = xc.open_dataset(
        path,
        decode_times=True,
        add_bounds=["X", "Y", "Z"],
    )
    try:
        return {
            "dims": dict(ds.sizes),
            "data_vars": sorted(str(key) for key in ds.data_vars.keys())[:5],
        }
    finally:
        ds.close()


def _open_with_xr_dataset(path: Path) -> dict[str, Any]:
    ds = xr.open_dataset(path, decode_times=True)
    try:
        return {
            "dims": dict(ds.sizes),
            "data_vars": sorted(str(key) for key in ds.data_vars.keys())[:5],
        }
    finally:
        ds.close()


OPENERS = {
    "xc.open_mfdataset": _open_with_xc_mfdataset,
    "xc.open_dataset": _open_with_xc_dataset,
    "xr.open_dataset": _open_with_xr_dataset,
}


def _worker(opener_name: str, path: str, iterations: int, result_queue: mp.Queue) -> None:
    opener = OPENERS[opener_name]
    file_path = Path(path)
    timings: list[float] = []
    sample: dict[str, Any] | None = None

    try:
        for _ in range(iterations):
            start = time.perf_counter()
            sample = opener(file_path)
            timings.append(time.perf_counter() - start)

        result_queue.put(
            {
                "status": "ok",
                "pid": os.getpid(),
                "iterations": iterations,
                "timings": timings,
                "sample": sample,
            }
        )
    except Exception as exc:  # pragma: no cover
        result_queue.put(
            {
                "status": "error",
                "pid": os.getpid(),
                "iterations": len(timings),
                "timings": timings,
                "error": repr(exc),
            }
        )


def _run_case(opener_name: str, path: Path, workers: int, iterations: int, timeout: int) -> dict[str, Any]:
    ctx = mp.get_context("fork")
    result_queue: mp.Queue = ctx.Queue()
    processes: list[mp.Process] = []
    started = time.perf_counter()

    for index in range(workers):
        proc = ctx.Process(
            target=_worker,
            args=(opener_name, str(path), iterations, result_queue),
            name=f"{opener_name}-worker-{index}",
        )
        proc.start()
        processes.append(proc)

    timed_out: list[int] = []
    exit_codes: dict[int, int | None] = {}

    for proc in processes:
        proc.join(timeout)
        if proc.is_alive():
            timed_out.append(proc.pid or -1)
            proc.terminate()
            proc.join(5)
            if proc.is_alive():
                proc.kill()
                proc.join(5)
        exit_codes[proc.pid or -1] = proc.exitcode

    results: list[dict[str, Any]] = []
    while True:
        try:
            results.append(result_queue.get_nowait())
        except queue.Empty:
            break

    ok_results = [item for item in results if item.get("status") == "ok"]
    error_results = [item for item in results if item.get("status") == "error"]
    all_timings = [timing for item in ok_results for timing in item.get("timings", [])]

    return {
        "opener": opener_name,
        "path": str(path),
        "realpath": str(path.resolve()),
        "workers": workers,
        "iterations_per_worker": iterations,
        "timeout_seconds": timeout,
        "elapsed_seconds": round(time.perf_counter() - started, 3),
        "completed_workers": len(ok_results),
        "errored_workers": len(error_results),
        "timed_out_workers": timed_out,
        "exit_codes": exit_codes,
        "total_successful_opens": len(all_timings),
        "min_open_seconds": round(min(all_timings), 6) if all_timings else None,
        "median_open_seconds": round(sorted(all_timings)[len(all_timings) // 2], 6)
        if all_timings
        else None,
        "max_open_seconds": round(max(all_timings), 6) if all_timings else None,
        "sample": ok_results[0].get("sample") if ok_results else None,
        "errors": error_results,
    }


def main() -> int:
    if not ANN_CLIMO_FILE.is_file():
        raise FileNotFoundError(f"Focused repro file not found: {ANN_CLIMO_FILE}")

    print(f"Focused repro file: {ANN_CLIMO_FILE}")
    print(f"Resolved source file: {ANN_CLIMO_FILE.resolve()}")
    print(
        f"Workers={NUM_WORKERS} Iterations/worker={NUM_ITERATIONS} Timeout={WORKER_TIMEOUT_SECONDS}s"
    )

    cases = [
        _run_case(
            opener_name=opener_name,
            path=ANN_CLIMO_FILE,
            workers=NUM_WORKERS,
            iterations=NUM_ITERATIONS,
            timeout=WORKER_TIMEOUT_SECONDS,
        )
        for opener_name in OPENERS
    ]

    print(json.dumps(cases, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())