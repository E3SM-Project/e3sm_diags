# Python 3.14 Runtime Risk Analysis (Consolidated)

This document consolidates prior prompt-response iterations into a single readable report.

## 1) Context and assumptions

- Repository: `e3sm_diags`
- Comparison: Python 3.13 (stable) vs Python 3.14 (intermittent stalls/cancels)
- Observed behavior: failures can occur at seemingly random variables/points, including with multiprocessing disabled
- Environment delta: mostly Python version and NumPy (`2.2.6` vs `2.4.2`)

## 2) Required scans executed

Commands used:

```bash
rg -n "multiprocessing|ProcessPoolExecutor|concurrent\.futures|dask|distributed|LocalCluster|Client|Nanny|WorkerPlugin|**del**|weakref\.finalize|atexit|faulthandler|signal|subprocess|resource|getrlimit|setrlimit|open_dataset|netCDF4|Dataset\(|h5py|ESMF|esmpy|xesmf|tempfile|NamedTemporaryFile|mkdtemp|cache|joblib|thread|ThreadPoolExecutor|asyncio" .
rg -n "PYTHON|OMP\_|OPENBLAS|MKL|NUMEXPR|HDF5|NETCDF|PROJ|GDAL|DASK|ESMF|MPI" .
rg -n "if **name** == ['\"]**main**['\"]|console_scripts|entry_points|argparse|click" .
```

Key outputs:

- Scan 1: `1284` matches (many in docs/notebooks; core runtime hotspots extracted below)
- Scan 2: env-var writes centered in `e3sm_diags/__init__.py:9-13`
- Scan 3: CLI entrypoints in `pyproject.toml:69-71` and `e3sm_diags/e3sm_diags_driver.py:552`

## 3) Executive summary (highest likelihood first)

1. Dataset/file lifecycle cleanup is not always deterministic in core paths.
2. Repeated `open_mfdataset()` + eager `.load(scheduler="sync")` patterns can amplify FD/memory pressure in loops.
3. Multiprocessing path has a brittle CPython-private bootstrap check (`_inheriting`).
4. ESMF/xESMF concurrency constraints are documented but not fully enforced by runtime knobs.
5. Exception swallowing can make runs appear to stop randomly.
6. Built-in fault/memory instrumentation is limited, which slows root-cause diagnosis.

## 4) Prioritized findings

### F1. GC-sensitive dataset cleanup in long-running flows

- Locations (examples):
  - `e3sm_diags/driver/utils/dataset_xr.py:423, 655, 1187, 1658`
  - `e3sm_diags/driver/arm_diags_driver.py:108`
  - `e3sm_diags/driver/precip_pdf_driver.py:124`
  - `e3sm_diags/driver/tc_analysis_driver.py:78, 441`
  - `e3sm_diags/driver/area_mean_time_series_driver.py:137`
- Why 3.14-sensitive: if cleanup timing shifts with GC/finalization behavior, native handles (HDF5/netCDF) can live longer than expected.
- Symptom fit:
  - Serial fail quickly: medium
  - Random variable stop points: high

### F2. Heavy I/O loops with `open_mfdataset()` and eager load

- Locations (examples):
  - `e3sm_diags/driver/utils/dataset_xr.py:1160-1197, 1490-1493, 1663-1709`
  - `e3sm_diags/driver/arm_diags_driver.py:82-153, 261-326`
  - `e3sm_diags/driver/precip_pdf_driver.py:176-320`
- Why 3.14-sensitive: repeated open/load/close pressure can surface latent lifecycle issues and memory spikes with slightly different runtime timing.
- Symptom fit:
  - Serial fail quickly: medium
  - Random variable stop points: high

### F3. CPython-private multiprocessing guard

- Location: `e3sm_diags/e3sm_diags_driver.py:405-414`
- Related scheduler config: `e3sm_diags/e3sm_diags_driver.py:357-369`
- Why 3.14-sensitive: private internals (like `_inheriting`) are not API-stable across CPython releases.
- Symptom fit:
  - Serial fail quickly: high when multiprocessing is enabled
  - Random variable stop points: medium

### F4. ESMF/xESMF runtime policy not fully enforced

- Documentation reference: `AGENTS.md:88-91`
- Concern: expected constraints (`processes=True`, `threads_per_worker=1`, ESMF isolation) are not globally enforced by a central runtime policy.
- Symptom fit:
  - Serial fail quickly: low-medium
  - Random variable stop points: medium-high

### F5. Exception swallowing can hide real failure boundaries

- Locations:
  - `e3sm_diags/parameter/core_parameter.py:407-413`
  - `e3sm_diags/run.py:104-109`
- Symptom fit:
  - Serial fail quickly: medium
  - Random variable stop points: medium

### F6. Missing first-class crash/memory instrumentation

- Observation: no central, always-available debug toggles around faulthandler/tracemalloc/resource snapshots.
- Impact: not likely causal, but strongly impedes diagnosis.

## 5) Finding matrix

| Finding | Primary location(s) | Risk mechanism | Suggested fix |
| --- | --- | --- | --- |
| GC-dependent xarray/netCDF cleanup | `dataset_xr.py`, `arm_diags_driver.py`, `tc_analysis_driver.py`, `area_mean_time_series_driver.py`, `precip_pdf_driver.py` | Native file handles/locks survive longer than intended | Explicit `close()` in `finally` or context managers |
| `open_mfdataset()` in loops | `dataset_xr.py`, `arm_diags_driver.py`, `precip_pdf_driver.py` | Repeated open/load cycles grow backend state and memory pressure | Open-load-close helpers; reuse safely when possible |
| Private multiprocessing guard | `e3sm_diags_driver.py:405-414` | Version-fragile use of private CPython internals | Replace with public API checks and explicit process-role logs |
| Dask orchestration observability gaps | `e3sm_diags_driver.py:357-369` | Worker creation/task failures lack robust diagnostics | Add task boundary logs, startup/timeout diagnostics |
| ESMF policy not centrally enforced | runtime + docs | Global-state native libs + unintended concurrency | Centralize process/thread policy for ESMF-sensitive paths |
| Exception swallowing | `core_parameter.py`, `run.py` | Failures appear as random stops | Add strict/debug mode and richer context in raised errors |
| Env vars only partially constrain thread pools | `e3sm_diags/__init__.py:9-13` | Oversubscription still possible (numexpr/blis/veclib) | Add centralized debug env policy |

## 6) Instrumentation patch plan (concise)

1. Process-start fault dumps
- File: `e3sm_diags/e3sm_diags_driver.py` near multiprocessing setup
- Add:
  - `faulthandler.enable(all_threads=True)`
  - `faulthandler.dump_traceback_later(600, repeat=True)`
  - `faulthandler.register(signal.SIGUSR1, all_threads=True)`

2. Memory/resource tracking
- File: `e3sm_diags/e3sm_diags_driver.py` near runtime bootstrap
- Add:
  - `tracemalloc.start(25)`
  - periodic RSS + open-FD logging (`resource.getrusage`, `/proc/self/fd` on Linux)
- File: `e3sm_diags/parameter/core_parameter.py` around set loop
- Add pre/post set snapshots

3. Structured dataset boundary logs
- File: `e3sm_diags/driver/utils/dataset_xr.py`
  - `_open_climo_dataset` (open start/end)
  - `_get_time_series_dataset_obj` (open/load/close)
  - `_subset_vars_and_load` (load/close)
- File: `e3sm_diags/driver/precip_pdf_driver.py`
  - `load_cached_pdf` (cache open/load/close)

4. Task boundary logs
- File: `e3sm_diags/e3sm_diags_driver.py`
- Add pre/post logs around serial and parallel task submission/execution boundaries

## 7) Minimal repro snippets (py313 vs py314)

### Repro 1: existing ESMF + Dask process script

```bash
python auxiliary_tools/cdat_regression_testing/933-esmf-mpi/mvce_grid_only.py
```

### Repro 2: core time-series open/load pressure

```python
import gc
import os

from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.parameter.core_parameter import CoreParameter

p = CoreParameter()
p.test_data_path = "/path/to/monthly_ts"
p.test_timeseries_input = True
p.test_start_yr = "1985"
p.test_end_yr = "2014"
p.sets = ["lat_lon"]
d = Dataset(p, data_type="test")

for i in range(200):
    ds = d.get_time_series_dataset("PRECT")
    _ = float(ds["PRECT"].mean().values)
    print(i, "fds", len(os.listdir("/proc/self/fd")))
    gc.collect()
```

### Repro 3: cache-open leak pattern

```python
import os
import xarray as xr

path = "/path/to/PRECT_PDF_global_test_..._ANN.nc"
for i in range(5000):
    ds = xr.open_dataset(path)
    _ = ds.attrs.get("start_year")
    if i % 100 == 0:
        print(i, len(os.listdir("/proc/self/fd")))
    # intentionally no close
```

Then repeat with `ds.load(); ds.close()` and compare stability.

## 8) Open-site audit (core runtime)

| Site | Location | Context manager | Explicit close | Escapes scope | In loops | Needs fix |
| --- | --- | --- | --- | --- | --- | --- |
| `xc.open_dataset` (time-slice full dataset) | `dataset_xr.py:484` | No | Yes (caller close at `dataset_xr.py:460`) | Temporary | Sometimes | No (patched) |
| `xc.open_mfdataset` (climo) | `dataset_xr.py:720,726` | No | Yes (`dataset_xr.py:670,407,1798`) | Yes | Yes | No (patched) |
| `xc.open_mfdataset` (time series) | `dataset_xr.py:1254` | No | Yes (`dataset_xr.py:1270,1798`) | Yes | Yes | No (patched) |
| `xr.open_dataset` (default mask in `Dataset`) | `dataset_xr.py:1734` | No | Yes (`dataset_xr.py:1735,1744`) | Returns loaded data | Sometimes | No (patched) |
| xarray helper wrapper | `tc_analysis_driver.py:41` | Yes | Yes | Returns loaded `DataArray` | Sometimes | No (patched) |
| direct xarray in TC path | `tc_analysis_driver.py:95,132,158,228` | Via helper | Yes | No | Sometimes | No (patched) |
| `netCDF4.Dataset` in obs metrics | `tc_analysis_driver.py:462` | Yes | Yes | No | Yes | No (patched) |
| `xr.open_dataset` cache load | `precip_pdf_driver.py:124` | No | No | Yes | Yes | Yes |
| `xr.open_dataset` mask load | `area_mean_time_series_driver.py:137` | No | No | Yes | Not inner loop | Yes |
| `xr.open_dataset` native full dataset | `dataset_native.py:171` | No | No | Yes | Potentially | Yes |
| `xr.open_dataset` ARM refs | `arm_diags_driver.py:108,198,278,357,418` | No | No | Local/repeated | Yes | Yes |

## 9) Why this can show up as Slurm cancel / BrokenProcessPool

- Open-FD growth and backend handle accumulation can stall I/O or hit system limits.
- Repeated open/load loops can increase peak RSS and trigger scheduler OOM/cancel behavior.
- Lock/handle lifetime variability from GC/finalization timing changes can produce non-deterministic stop points.
- In parallel paths, worker death from the same resource pressure can surface as `BrokenProcessPool` or cancelled futures.

## 10) Minimal patch set summary (already drafted)

High-impact patch themes already proposed:

- Add deterministic close wrappers in `tc_analysis_driver.py` for xarray/netCDF access.
- Add `finally`-based close in key `dataset_xr.py` open/load paths.
- Add optional env-gated resource instrumentation (`E3SM_DIAGS_DEBUG_RESOURCES=1`) around dataset boundaries.

If needed, keep the full git-style diff in a separate file (`codex-analysis-diff.md`) so this report remains concise.

## 11) Validation note

Full tests were not executed in the shell used for this analysis because the available interpreter was incompatible with this repo configuration. Validate in the project conda environment.

## Appendix A: Original prompts (for provenance)

<details>
<summary>Show original prompt text</summary>

The original prompt chain requested:

- repository-wide scans for multiprocessing/dask/cleanup/runtime env vars,
- a Python 3.14 runtime-sensitivity analysis,
- a risk-prioritized findings list,
- concrete low-risk mitigations,
- instrumentation insertion points,
- and minimal repro scripts.

A follow-up prompt requested:

- full open-site audit table,
- rationale mapping to `Slurm cancelled` / `BrokenProcessPool`,
- and a unified diff for minimal cleanup/instrumentation patches.

</details>
