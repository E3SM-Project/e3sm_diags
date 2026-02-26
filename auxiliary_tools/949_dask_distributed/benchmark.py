"""Benchmarking harness for comparing Dask scheduler configurations.

This script benchmarks E3SM diagnostics execution across different Dask
scheduler configurations:
  1. dask.bag with "processes" scheduler (legacy, default)
  2. dask.distributed with LocalCluster

It records wall time and peak memory usage (via tracemalloc), then
generates comparison bar charts using matplotlib.

Usage
-----
This script is designed to be run on an HPC system where the required
data paths are available. It must be executed as:

    python -m auxiliary_tools.949_dask_distributed.benchmark \\
        --results-dir /path/to/results \\
        --test-data-path /path/to/test/data \\
        --reference-data-path /path/to/ref/data \\
        --num-workers 4

All arguments are optional and have defaults pointing to common LCRC paths.

Notes
-----
- This script does NOT modify production execution paths.
- Results are written to ``<results_dir>/benchmark_results/``.
- The benchmark is isolated and safe to run without affecting other code.
"""

import argparse
import json
import os
import sys
import time
import tracemalloc
from typing import Any

# Attempt to import matplotlib; if unavailable, plotting is skipped.
try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


def _create_parameter(args: argparse.Namespace) -> Any:
    """Create a CoreParameter configured for benchmarking.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    CoreParameter
        A configured parameter object.
    """
    from e3sm_diags.parameter.core_parameter import CoreParameter

    param = CoreParameter()
    param.test_data_path = args.test_data_path
    param.reference_data_path = args.reference_data_path
    param.test_name = args.test_name
    param.run_type = "model_vs_obs"
    param.diff_title = "Model - Observations"
    param.output_format = ["png"]
    param.output_format_subplot = []
    param.multiprocessing = True
    param.num_workers = args.num_workers
    param.save_netcdf = False

    return param


def _run_benchmark_bag(
    args: argparse.Namespace,
) -> dict[str, float]:
    """Run diagnostics using the dask.bag 'processes' scheduler.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict[str, float]
        Dictionary with 'wall_time_s' and 'peak_memory_mb'.
    """
    from e3sm_diags.run import runner

    param = _create_parameter(args)
    param.results_dir = os.path.join(args.results_dir, "benchmark_bag")
    param.dask_scheduler_type = "processes"

    runner.sets_to_run = args.sets.split(",")

    tracemalloc.start()
    t0 = time.perf_counter()

    try:
        runner.run_diags([param], use_cfg=False)
    except Exception as e:
        print(f"Bag benchmark failed: {e}")

    wall_time = time.perf_counter() - t0
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "wall_time_s": round(wall_time, 2),
        "peak_memory_mb": round(peak_mem / (1024 * 1024), 2),
    }


def _run_benchmark_distributed(
    args: argparse.Namespace,
) -> dict[str, float]:
    """Run diagnostics using the dask.distributed scheduler.

    Parameters
    ----------
    args : argparse.Namespace
        Parsed command-line arguments.

    Returns
    -------
    dict[str, float]
        Dictionary with 'wall_time_s' and 'peak_memory_mb'.
    """
    from dask.distributed import LocalCluster

    from e3sm_diags.run import runner

    param = _create_parameter(args)
    param.results_dir = os.path.join(args.results_dir, "benchmark_distributed")
    param.dask_scheduler_type = "distributed"
    param.dask_memory_limit = args.memory_limit

    runner.sets_to_run = args.sets.split(",")

    cluster = LocalCluster(
        n_workers=args.num_workers,
        threads_per_worker=1,
        processes=True,
        memory_limit=args.memory_limit,
        resources={"ESMF": 1},
    )

    tracemalloc.start()
    t0 = time.perf_counter()

    try:
        runner.run_diags([param], use_cfg=False, dask_cluster=cluster)
    except Exception as e:
        print(f"Distributed benchmark failed: {e}")

    wall_time = time.perf_counter() - t0
    _, peak_mem = tracemalloc.get_traced_memory()
    tracemalloc.stop()

    return {
        "wall_time_s": round(wall_time, 2),
        "peak_memory_mb": round(peak_mem / (1024 * 1024), 2),
    }


def _generate_plots(
    results: dict[str, dict[str, float]], output_dir: str
) -> None:
    """Generate comparison bar charts from benchmark results.

    Parameters
    ----------
    results : dict[str, dict[str, float]]
        Benchmark results keyed by scheduler name.
    output_dir : str
        Directory to write plot files.
    """
    if not HAS_MATPLOTLIB:
        print("matplotlib not available; skipping plot generation.")
        return

    schedulers = list(results.keys())
    wall_times = [results[s]["wall_time_s"] for s in schedulers]
    peak_mems = [results[s]["peak_memory_mb"] for s in schedulers]

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))

    # Wall time comparison
    axes[0].bar(schedulers, wall_times, color=["steelblue", "coral"])
    axes[0].set_ylabel("Wall Time (seconds)")
    axes[0].set_title("Wall Time Comparison")
    for i, v in enumerate(wall_times):
        axes[0].text(i, v + 0.5, f"{v:.1f}s", ha="center", fontsize=9)

    # Peak memory comparison
    axes[1].bar(schedulers, peak_mems, color=["steelblue", "coral"])
    axes[1].set_ylabel("Peak Memory (MB)")
    axes[1].set_title("Peak Memory Comparison")
    for i, v in enumerate(peak_mems):
        axes[1].text(i, v + 0.5, f"{v:.1f}MB", ha="center", fontsize=9)

    plt.tight_layout()
    plot_path = os.path.join(output_dir, "benchmark_comparison.png")
    plt.savefig(plot_path, dpi=150)
    plt.close()
    print(f"Benchmark plot saved to: {plot_path}")


def main() -> None:
    """Entry point for the benchmarking harness."""
    parser = argparse.ArgumentParser(
        description="Benchmark E3SM diagnostics with different Dask schedulers."
    )
    parser.add_argument(
        "--results-dir",
        default="benchmark_results",
        help="Base directory for benchmark outputs.",
    )
    parser.add_argument(
        "--test-data-path",
        default="",
        help="Path to test (model) data.",
    )
    parser.add_argument(
        "--reference-data-path",
        default="",
        help="Path to reference (obs) data.",
    )
    parser.add_argument(
        "--test-name",
        default="benchmark_test",
        help="Name for the test model run.",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=4,
        help="Number of Dask workers.",
    )
    parser.add_argument(
        "--memory-limit",
        default="auto",
        help="Memory limit per worker (e.g., '2GB', 'auto').",
    )
    parser.add_argument(
        "--sets",
        default="lat_lon",
        help="Comma-separated list of diagnostic sets to run.",
    )
    parser.add_argument(
        "--skip-bag",
        action="store_true",
        help="Skip the dask.bag benchmark.",
    )
    parser.add_argument(
        "--skip-distributed",
        action="store_true",
        help="Skip the dask.distributed benchmark.",
    )

    args = parser.parse_args()

    output_dir = os.path.join(args.results_dir, "benchmark_results")
    os.makedirs(output_dir, exist_ok=True)

    results: dict[str, dict[str, float]] = {}

    if not args.skip_bag:
        print("\n" + "=" * 60)
        print("Running benchmark: dask.bag (processes)")
        print("=" * 60)
        results["bag"] = _run_benchmark_bag(args)
        print(f"  Wall time: {results['bag']['wall_time_s']}s")
        print(f"  Peak memory: {results['bag']['peak_memory_mb']}MB")

    if not args.skip_distributed:
        print("\n" + "=" * 60)
        print("Running benchmark: dask.distributed")
        print("=" * 60)
        results["distributed"] = _run_benchmark_distributed(args)
        print(f"  Wall time: {results['distributed']['wall_time_s']}s")
        print(
            f"  Peak memory: {results['distributed']['peak_memory_mb']}MB"
        )

    # Save results as JSON
    results_path = os.path.join(output_dir, "benchmark_results.json")
    with open(results_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nResults saved to: {results_path}")

    # Generate comparison plots
    if len(results) >= 2:
        _generate_plots(results, output_dir)

    # Print summary
    print("\n" + "=" * 60)
    print("Benchmark Summary")
    print("=" * 60)
    for name, metrics in results.items():
        wall = metrics['wall_time_s']
        mem = metrics['peak_memory_mb']
        print(f"  {name}: {wall}s, {mem}MB peak")

    if "bag" in results and "distributed" in results:
        speedup = results["bag"]["wall_time_s"] / max(
            results["distributed"]["wall_time_s"], 0.01
        )
        print(f"\n  Distributed vs Bag speedup: {speedup:.2f}x")


if __name__ == "__main__":
    main()
