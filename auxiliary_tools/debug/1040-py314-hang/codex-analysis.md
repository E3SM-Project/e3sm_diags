## Prompt for Codex
```markdown
Run these commands first (capture output in your analysis):

1) Repository-wide scans:
rg -n "multiprocessing|ProcessPoolExecutor|concurrent\.futures|dask|distributed|LocalCluster|Client|Nanny|WorkerPlugin|__del__|weakref\.finalize|atexit|faulthandler|signal|subprocess|resource|getrlimit|setrlimit|open_dataset|netCDF4|Dataset\(|h5py|ESMF|esmpy|xesmf|tempfile|NamedTemporaryFile|mkdtemp|cache|joblib|thread|ThreadPoolExecutor|asyncio" .

2) Where runtime env vars are referenced:
rg -n "PYTHON|OMP_|OPENBLAS|MKL|NUMEXPR|HDF5|NETCDF|PROJ|GDAL|DASK|ESMF|MPI" .

3) If there is a driver/CLI entrypoint, locate it:
rg -n "if __name__ == ['\"]__main__['\"]|console_scripts|entry_points|argparse|click" .

Now analyze the e3sm_diags codebase for Python 3.14 runtime-related failure risks (hangs, random stops, BrokenProcessPool, unexpected termination), focusing on code paths sensitive to CPython runtime changes rather than pure package-version issues.

Context to assume:
- Python 3.13 runs complete; Python 3.14 sometimes stalls/cancels even with multiprocessing disabled.
- Failures occur at seemingly random points/variables.
- Conda environments differ mainly by Python 3.13 vs 3.14 and NumPy 2.2.6 vs 2.4.2; most other package versions are the same.

Tasks:
1) Identify all places in the repo that interact with:
   - multiprocessing / concurrent.futures (ProcessPoolExecutor), start methods, pools, Manager, shared_memory
   - dask/distributed process orchestration (LocalCluster, Client, worker lifecycle, nanny, plugins)
   - thread/process signals, timeouts, subprocess, os.fork, resource limits, atexit handlers
   - object finalization patterns: __del__, weakref.finalize, reliance on GC for cleanup
   - context management and explicit resource cleanup (closing datasets/files, temp dirs, caches)
   - heavy I/O loops that repeatedly open/close datasets (xarray/netCDF4), temp files, caching directories
   - environment variables that affect runtime (PYTHON*, OMP_*, OPENBLAS_*, MKL_*, HDF5_*, PROJ_*, GDAL_*, DASK_*, ESMF_*, MPI_*)
2) For each finding, report:
   - file path + line numbers
   - what runtime-sensitive behavior exists and why it might differ under Python 3.14 (GC/finalization timing, signal handling, spawn semantics, interpreter shutdown behavior, etc.)
   - whether it can plausibly explain (a) “serial run fails quickly” and/or (b) “random variable stop points”
3) Propose concrete mitigations, preferring low-risk changes:
   - enforce explicit close()/cleanup; convert to context-manager usage; ensure dask clients/clusters are closed
   - avoid relying on __del__ / GC for releasing native resources
   - set/centralize multiprocessing start method selection (if used) and ensure it’s set early
   - reduce thread oversubscription (OMP/BLAS) if plausible; note where to set env vars
   - add robust timeouts and better error/trace capture around worker creation and I/O loops
4) Provide a short “instrumentation patch plan”:
   - where to add faulthandler.enable(), faulthandler.dump_traceback_later(), signal-triggered dumps
   - where to add tracemalloc and periodic RSS logging
   - where to add structured logging around dataset open/load/close and task boundaries
5) If possible, derive 1–3 minimal repro snippets from the repo (small scripts that exercise the same code paths) that could be run under py313 vs py314.

Output format:
- A prioritized list of findings (highest likelihood first)
- A table: Finding | Location | Risk mechanism | Suggested fix
- A concise instrumentation plan (files + exact insertion points)
- Minimal repro snippets (if feasible)
```

## Overview

This analysis focuses on Python 3.14 runtime-sensitivity risks in `e3sm_diags`
that can explain non-deterministic hangs, cancellation, or random stop points
even when multiprocessing is disabled.

Highest-probability contributors are:

1. Dataset/file lifecycles that rely on GC timing instead of deterministic close.
2. Repeated heavy I/O loops (`open_mfdataset` + eager `load`) without explicit
   cleanup boundaries.
3. Parallel runtime paths with brittle process-bootstrap assumptions and limited
   worker/task-level diagnostics.

The report below includes:

1. Verified repository scans and exact code locations.
2. Prioritized findings with py3.14-sensitive risk mechanisms.
3. Low-risk mitigations emphasizing explicit cleanup and safer runtime control.
4. A focused instrumentation plan and minimal repro snippets for py313 vs py314.

**Requested scans executed**
I ran your three exact `rg` commands first. Key outputs:

- Command 1 returned `1284` matches (many in docs/notebooks; core runtime hotspots identified below).
- Command 2 found runtime env-var writes primarily in package init: [**init**.py:9-13](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/__init__.py:9).
- Command 3 found CLI entrypoints at [pyproject.toml:69-71](/home/ac.tvo/E3SM-Project/e3sm_diags/pyproject.toml:69) and [e3sm_diags_driver.py:552](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:552).

**Prioritized findings (highest likelihood first)**

1. **Core runtime relies on GC/finalization for dataset/file cleanup (no explicit `close()` in core path).**  
   Locations: dataset open sites in [dataset_xr.py:423](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:423), [dataset_xr.py:655](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:655), [dataset_xr.py:1187](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1187), [dataset_xr.py:1658](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1658), [arm_diags_driver.py:108](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/arm_diags_driver.py:108), [precip_pdf_driver.py:124](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/precip_pdf_driver.py:124), [tc_analysis_driver.py:78](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py:78), [area_mean_time_series_driver.py:137](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/area_mean_time_series_driver.py:137).  
   Also `netCDF4.Dataset` opened without context manager in [tc_analysis_driver.py:441](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py:441).  
   Plausibility for your symptoms: **(a) medium, (b) high**.  
   Inference: under Python 3.14, changed GC/finalization timing can expose latent file-handle/lock lifetime bugs that looked harmless before.

2. **Heavy repeated I/O loops with `open_mfdataset()` + eager `.load(scheduler="sync")` but no explicit source close.**  
   Locations: [dataset_xr.py:1160-1197](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1160), [dataset_xr.py:1663-1709](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1663), [dataset_xr.py:1490-1493](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1490), and multi-loop drivers like [arm_diags_driver.py:82-153](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/arm_diags_driver.py:82), [arm_diags_driver.py:261-326](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/arm_diags_driver.py:261), [precip_pdf_driver.py:176-320](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/precip_pdf_driver.py:176).  
   Plausibility: **(a) medium, (b) high**.

3. **Multiprocessing path uses brittle CPython-private bootstrap check (`_inheriting`).**  
   Locations: [e3sm_diags_driver.py:405-414](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:405), with Dask process scheduler config in [e3sm_diags_driver.py:357-369](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:357).  
   Plausibility: **(a) high when multiprocessing on, (b) medium**.  
   Inference: private multiprocessing internals can change across CPython versions; this check is version-fragile.

4. **ESMF/xESMF thread/process constraints are documented but not enforced in runtime config.**  
   Documentation: [AGENTS.md:88-91](/home/ac.tvo/E3SM-Project/e3sm_diags/AGENTS.md:88).  
   Core sets only Dask `processes` + `forkserver` locally, but no runtime enforcement for `threads_per_worker=1` or ESMF resource isolation.  
   Plausibility: **(a) low-medium, (b) medium-high** for random stop points in regridding-heavy sets.

5. **Exception swallowing can make failures look random/incomplete instead of hard-failing with clear traceback boundaries.**  
   Locations: [core_parameter.py:407-413](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py:407), [run.py:104-109](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/run.py:104).  
   Plausibility: **(a) medium, (b) medium** for “random stop points” perception.

6. **No built-in fault/memory instrumentation (`faulthandler`, signal-triggered dumps, tracemalloc, RSS).**  
   `rg` found no active usage in core runtime.  
   Plausibility: not causal, but strongly blocks diagnosis.

**Finding table**

| Finding                                           | Location                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                | Risk mechanism                                                                              | Suggested fix                                                                                                                |
| ------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | -------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------- |
| GC-dependent cleanup of xarray/netCDF handles     | [dataset_xr.py:423](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:423), [dataset_xr.py:655](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:655), [dataset_xr.py:1187](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1187), [arm_diags_driver.py:108](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/arm_diags_driver.py:108), [tc_analysis_driver.py:78](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py:78) | Delayed finalization can retain HDF5/netCDF file locks/fds, causing non-deterministic stalls | Use context managers or explicit `close()` in `finally`; after `load()`, detach/copy needed arrays and close source datasets |
| `open_mfdataset()` in repeated loops              | [dataset_xr.py:1160](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1160), [dataset_xr.py:1663](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1663), [precip_pdf_driver.py:176](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/precip_pdf_driver.py:176)                                                                                                                                                                                                               | Repeated lazy-open + sync load without deterministic close can accumulate backend state      | Add helper to open/load/close atomically; reuse datasets where possible per var/season                                       |
| Private multiprocessing bootstrap guard           | [e3sm_diags_driver.py:413](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:413)                                                                                                                                                                                                                                                                                                                                                                                                                                    | CPython internal attr `_inheriting` is not API-stable; may break in 3.14+                    | Replace with safer guard strategy; avoid private attr reliance; add explicit process role logging                            |
| Dask process orchestration lacks robustness knobs | [e3sm_diags_driver.py:357-369](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:357)                                                                                                                                                                                                                                                                                                                                                                                                                                | Worker startup/task failures can surface as pool breaks/cancel without rich context          | Centralize start-method policy, add startup/task boundary logs, timeout wrappers, and worker error dumps                     |
| ESMF/xESMF safety constraints not enforced        | [AGENTS.md:88-91](/home/ac.tvo/E3SM-Project/e3sm_diags/AGENTS.md:88)                                                                                                                                                                                                                                                                                                                                                                                                                                                                    | Global-state native libs can deadlock/crash with unintended concurrency                      | Enforce `threads_per_worker=1` equivalent in all parallel paths; keep ESMF work process-isolated                             |
| `netCDF4.Dataset` unclosed                        | [tc_analysis_driver.py:441](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/tc_analysis_driver.py:441)                                                                                                                                                                                                                                                                                                                                                                                                                           | Native resource leak across looped file opens                                                | `with netcdffile(...) as nc:`                                                                                                |
| Cache file opens without close                    | [precip_pdf_driver.py:124](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/precip_pdf_driver.py:124)                                                                                                                                                                                                                                                                                                                                                                                                                             | Cache hit path can leak file descriptors over many seasons/refs                              | `xr.open_dataset(...); ds.load(); ds.close(); return loaded copy`                                                            |
| Exception swallowing obscures root cause          | [core_parameter.py:407](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py:407), [run.py:104](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/run.py:104)                                                                                                                                                                                                                                                                                                                                                      | Run appears to “stop randomly” when failures are suppressed                                  | Optional strict mode default in debug runs; include set/var context in raised errors                                         |
| Env vars only partially constrain threading       | [**init**.py:9-13](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/__init__.py:9)                                                                                                                                                                                                                                                                                                                                                                                                                                                       | Some libs (e.g. numexpr/blis/veclib) may still oversubscribe                                 | Set/centralize `NUMEXPR_NUM_THREADS`, `VECLIB_MAXIMUM_THREADS`, `BLIS_NUM_THREADS` for debug runs                            |

**Concise instrumentation patch plan**

1. **Enable crash dumps and signal-triggered tracebacks at process start**  
   Insert in [e3sm_diags_driver.py:405-416](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:405):

- `faulthandler.enable(all_threads=True)`
- `faulthandler.dump_traceback_later(600, repeat=True)`
- `faulthandler.register(signal.SIGUSR1, all_threads=True)`

2. **Add tracemalloc + RSS/open-fd periodic logging**  
   Insert near [e3sm_diags_driver.py:416](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:416):

- `tracemalloc.start(25)`
- periodic logger (every N sec) capturing RSS (`resource.getrusage`) + fd count (`/proc/self/fd` on Linux).
- log snapshots before/after each set in [core_parameter.py:388-406](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/parameter/core_parameter.py:388).

3. **Structured I/O boundary logs**  
   Add start/end timing + file path + var/set in:

- [dataset_xr.py:608](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:608) `_open_climo_dataset`
- [dataset_xr.py:1168](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1168) `_get_time_series_dataset_obj`
- [dataset_xr.py:1663](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/utils/dataset_xr.py:1663) `_subset_vars_and_load`
- [precip_pdf_driver.py:80](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/driver/precip_pdf_driver.py:80) `load_cached_pdf`

4. **Task boundary logs for parallel execution**  
   Add detailed pre/post logs around [e3sm_diags_driver.py:367-369](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:367) and around each serial run in [e3sm_diags_driver.py:322-324](/home/ac.tvo/E3SM-Project/e3sm_diags/e3sm_diags/e3sm_diags_driver.py:322).

**Minimal repro snippets (py313 vs py314)**

1. **Existing ESMF+Dask process repro (already in repo)**  
   Run:  
   `python auxiliary_tools/cdat_regression_testing/933-esmf-mpi/mvce_grid_only.py`  
   File: [mvce_grid_only.py](/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/cdat_regression_testing/933-esmf-mpi/mvce_grid_only.py)

2. **Dataset open/load/close stress (core code path)**

```python
import os, gc
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.driver.utils.dataset_xr import Dataset

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

Compare fd growth/hang behavior between 3.13 and 3.14.

3. **Cache-open leak repro (precip PDF path)**

```python
import xarray as xr, os
path = "/path/to/PRECT_PDF_global_test_..._ANN.nc"
for i in range(5000):
    ds = xr.open_dataset(path)
    _ = ds.attrs.get("start_year")
    if i % 100 == 0:
        print(i, len(os.listdir("/proc/self/fd")))
    # intentionally no close, mirrors risky pattern
```

Then rerun with `ds.load(); ds.close()` to confirm stabilization.

If you want, I can turn this into a concrete patch with explicit close/cleanup and env-gated instrumentation in the exact files above.
