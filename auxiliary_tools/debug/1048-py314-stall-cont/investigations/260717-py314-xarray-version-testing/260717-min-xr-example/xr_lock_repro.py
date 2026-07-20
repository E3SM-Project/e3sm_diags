#!/usr/bin/env python3
"""Stress Xarray's default NetCDF4 backend lock with a NetCDF3 file.

This script approximates the E3SM Diags climatology lifecycle without importing
E3SM Diags or xCDAT. Each forked worker keeps test, reference, and mask datasets
open concurrently, synchronously loads derived and mask data, deep-copies the
results, closes the files, and repeats the sequence.

The stall is intermittent. Run multiple rounds with Xarray 2026.1.0 or newer
under Python 3.14.3, then compare against ``--lock false``.

Examples
--------
Run with Xarray's default NetCDF4 backend lock::
    srun --pty --nodes=1 --time=01:00:00 /bin/bash
    conda activate ed_1048_xr_latest_2026070_py3143
    python auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/260717-min-xr-example/xr_lock_repro.py

Run the lock-free control::

    python auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/260717-min-xr-example/xr_lock_repro.py --lock false
"""
from __future__ import annotations

import argparse
import faulthandler
import multiprocessing as mp
import os
import signal
import sys
import time
from pathlib import Path
from typing import Any, Literal

import xarray as xr


DEFAULT_FILE = (
    Path(__file__).resolve().parents[2]
    / "archive"
    / "xarray-min"
    / "v3.LR.historical_0051_ANN_198501_201412_climo.nc"
)
_STACK_FILE: Any = None


def _initialize_worker(stack_dir: str) -> None:
    """Register a signal-triggered stack dump for a worker process."""
    global _STACK_FILE

    stack_path = Path(stack_dir) / f"worker-{os.getpid()}.log"
    _STACK_FILE = stack_path.open("w", encoding="utf-8")
    faulthandler.register(signal.SIGUSR1, file=_STACK_FILE, all_threads=True)


def _open_kwargs(lock_mode: Literal["default", "false"]) -> dict[str, Any]:
    kwargs: dict[str, Any] = {
        "engine": "netcdf4",
        "decode_times": True,
    }
    if lock_mode == "false":
        kwargs["lock"] = False

    return kwargs


def _load_albedo(ds: xr.Dataset) -> xr.Dataset:
    """Derive and detach ALBEDO, matching the relevant E3SM Diags read."""
    if "FSNTOA" not in ds or "SOLIN" not in ds:
        available = ", ".join(sorted(str(name) for name in ds.data_vars))
        raise KeyError(f"FSNTOA and SOLIN are required; found: {available}")

    result = (ds["FSNTOA"] / ds["SOLIN"]).to_dataset(name="ALBEDO")
    result.load(scheduler="sync")
    return result.copy(deep=True)


def _load_mask(ds: xr.Dataset) -> xr.Dataset:
    """Load and detach available land/ocean mask variables."""
    mask_vars = [name for name in ("LANDFRAC", "OCNFRAC") if name in ds]
    if not mask_vars:
        mask_vars = [str(next(iter(ds.data_vars)))]

    result = ds[mask_vars]
    result.load(scheduler="sync")
    return result.copy(deep=True)


def _worker_task(args: tuple[int, str, str, int]) -> dict[str, Any]:
    task_id, filepath, lock_mode, opens_per_task = args
    kwargs = _open_kwargs(lock_mode)  # type: ignore[arg-type]
    started = time.perf_counter()

    for _ in range(opens_per_task):
        ds_test: xr.Dataset | None = None
        ds_ref: xr.Dataset | None = None
        ds_mask: xr.Dataset | None = None
        try:
            # Keep all three backend file managers alive while data are loaded.
            ds_test = xr.open_mfdataset(
                paths=filepath,
                coords="minimal",
                compat="override",
                **kwargs,
            )
            ds_ref = xr.open_mfdataset(
                paths=filepath,
                coords="minimal",
                compat="override",
                **kwargs,
            )
            ds_mask = xr.open_dataset(filepath, **kwargs)

            _load_albedo(ds_test)
            _load_albedo(ds_ref)
            _load_mask(ds_mask)
        finally:
            for ds in (ds_mask, ds_ref, ds_test):
                if ds is not None:
                    ds.close()

    return {
        "task_id": task_id,
        "pid": os.getpid(),
        "elapsed_seconds": round(time.perf_counter() - started, 3),
    }


def _dump_worker_stacks(pool: Any) -> list[int]:
    pids: list[int] = []
    for process in pool._pool:
        if process.is_alive() and process.pid is not None:
            pids.append(process.pid)
            os.kill(process.pid, signal.SIGUSR1)

    return pids


def _run_round(
    *,
    round_number: int,
    filepath: Path,
    lock_mode: Literal["default", "false"],
    workers: int,
    tasks: int,
    opens_per_task: int,
    timeout: int,
    stack_dir: Path,
) -> bool:
    ctx = mp.get_context("fork")
    pool = ctx.Pool(
        processes=workers,
        initializer=_initialize_worker,
        initargs=(str(stack_dir),),
    )
    task_args = [
        (task_id, str(filepath), lock_mode, opens_per_task)
        for task_id in range(tasks)
    ]
    started = time.perf_counter()

    try:
        pending = pool.map_async(_worker_task, task_args, chunksize=1)
        results = pending.get(timeout=timeout)
    except mp.TimeoutError:
        stuck_pids = _dump_worker_stacks(pool)
        pool.terminate()
        pool.join()
        print(
            f"Round {round_number}: TIMEOUT after {timeout}s; "
            f"worker stacks requested for {stuck_pids}",
            flush=True,
        )
        return False
    except BaseException:
        pool.terminate()
        pool.join()
        raise

    pool.close()
    pool.join()
    elapsed = time.perf_counter() - started
    slowest = max(float(result["elapsed_seconds"]) for result in results)
    print(
        f"Round {round_number}: completed {len(results)}/{tasks} tasks in "
        f"{elapsed:.2f}s (slowest task: {slowest:.2f}s)",
        flush=True,
    )
    return True


def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--file", type=Path, default=DEFAULT_FILE)
    parser.add_argument("--lock", choices=("default", "false"), default="default")
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--rounds", type=int, default=10)
    parser.add_argument("--tasks", type=int, default=80)
    parser.add_argument("--opens-per-task", type=int, default=3)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument(
        "--stack-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "worker-stacks",
    )
    return parser.parse_args()


def main() -> int:
    args = _parse_args()
    filepath = args.file.resolve()
    stack_dir = args.stack_dir.resolve()

    if not filepath.is_file():
        raise FileNotFoundError(f"NetCDF file not found: {filepath}")
    if min(args.workers, args.rounds, args.tasks, args.opens_per_task, args.timeout) < 1:
        raise ValueError("Worker, round, task, open, and timeout values must be positive")

    stack_dir.mkdir(parents=True, exist_ok=True)
    print(f"Python: {sys.version.split()[0]}")
    print(f"Xarray: {xr.__version__}")
    print(f"File: {filepath}")
    print(f"Lock: {args.lock}")
    print(
        f"Workers={args.workers} rounds={args.rounds} tasks/round={args.tasks} "
        f"opens/task={args.opens_per_task} timeout={args.timeout}s"
    )
    print(f"Worker stacks: {stack_dir}")

    for round_number in range(1, args.rounds + 1):
        completed = _run_round(
            round_number=round_number,
            filepath=filepath,
            lock_mode=args.lock,
            workers=args.workers,
            tasks=args.tasks,
            opens_per_task=args.opens_per_task,
            timeout=args.timeout,
            stack_dir=stack_dir,
        )
        if not completed:
            return 124

    print(f"Completed all {args.rounds} rounds without a stall")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())