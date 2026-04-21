## Python 3.14 Stall Notes

## TL;DR

- I can still reproduce intermittent stalls on Python 3.14 with
  `xarray>=2026.01.0` using the NERSC minimum script with `num_workers=8`.
- The strongest evidence points to the xarray/netCDF4 backend lock-acquire
  path during climatology open/read work.
- The stall still happens after both of these changes:
  - switching single-file climatologies from `xc.open_mfdataset()` to
    `xc.open_dataset()`
  - switching multiprocessing from `fork` to `forkserver`
- A 10-run diagnostic with `lock=False` in the climo open path completed
  without stalling, while earlier `lock=True` repros typically hung by run 2.
  That is strong evidence that the backend lock path is required for
  reproduction.
- `lock=False` is evidence, not the target fix.

## What Is Happening

Some parallel runs finish and some hang.

When the hang occurs, the stuck worker is not spending time in plotting,
regridding, or later output steps. The best worker traceback so far points to
the xarray/netCDF4 backend lock path during climatology open/read work.

Representative stack at hang time:

```text
xarray/backends/locks.py:66  __enter__
xarray/backends/file_manager.py:181  _optional_lock
xarray/backends/file_manager.py:217  _acquire_with_cache_info
xarray/backends/netCDF4_.py:532  _acquire
xarray/backends/netCDF4_.py:108  get_array
```

## What We Tested

### 1. Targeted stuck-worker traceback

- A targeted `SIGUSR1` removed the ambiguity about which worker was actually
  stuck.
- In the first clear case:
  - stuck worker: `pid=1313623`
  - task context: test `ANN` `ALBEDO` climo file
    `v3.LR.historical_0051_ANN_198501_201412_climo.nc`
  - last `e3sm_diags` log line from that worker:
    `Climo backend open_mfdataset start`

Conclusion: the hang is tied to the climo backend open/read path.

### 2. Single-file `open_dataset()` instead of `open_mfdataset()`

I changed single concrete climatology paths to use `xc.open_dataset()` and kept
`xc.open_mfdataset()` only for true multi-file or globbed inputs.

Result:

- the stall still reproduced
- the later stuck worker showed the last log line
  `Climo backend open_dataset start`
- the worker traceback still showed the same xarray/netCDF4 lock path

Conclusion: this is not just a single-file misuse of `open_mfdataset()`.

### 3. `forkserver` instead of `fork`

I re-ran with `forkserver`.

Result:

- the stall still reproduced
- the worker shown finishing in the pasted excerpt (`pid=2254797`) was not the
  stuck worker
- the stuck worker traceback sink was `pid=2254799`
- that traceback still showed the same xarray/netCDF4 lock-acquire path

Conclusion: plain post-`fork` inherited state is not required for
reproduction.

### 4. `lock=False` as a diagnostic

I re-ran the same looped ALBEDO repro with `lock=False` passed through the
xarray/xCDAT open calls.

Result:

- all 10 `lock=False` runs completed without stalling
- earlier `lock=True` repros usually stalled by run 2

Conclusion: the xarray/netCDF4 backend lock path appears to be required for
reproduction.

### 5. Focused open-only repro

I ran a focused concurrent-open repro against the exact `ANN` climo file with
8 worker processes and 20 iterations per worker.

Result:

- `xc.open_mfdataset()`: 160/160 successful opens, no timeouts
- `xc.open_dataset()`: 160/160 successful opens, no timeouts
- `xr.open_dataset()`: 160/160 successful opens, no timeouts

Conclusion: simple concurrent open-only activity is not enough to trigger the
hang. The failure likely depends on more of the real `e3sm_diags` lifecycle.

### 6. `h5netcdf` engine check

I also tested whether `h5netcdf` could serve as a backend workaround for the
same staged climatology file.

Result:

- the open failed immediately with
  `OSError: Unable to synchronously open file (file signature not found)`
- `ncdump -k` on the same file reported `64-bit offset`

Conclusion: these staged climo files are not HDF5-backed netCDF4 files, so
`h5netcdf` is not currently a usable drop-in fallback for this dataset.

## What This Means Right Now

### Strongest current interpretation

- The hang is still strongly tied to the climo backend open path under
  multiprocessing.
- The stuck worker stack shows blocking in xarray/netCDF4 lock acquisition.
- The issue survives both of the obvious reductions:
  - `open_dataset()` for concrete single-file climatologies
  - `forkserver` instead of `fork`
- The issue does not show up in a simplified open-only repro.
- The obvious alternate-engine fallback is not currently usable on these
  staged climo files, because they are `64-bit offset` rather than
  HDF5-backed netCDF4.

### Practical implication

The trigger is probably not just "opening a file in parallel".

It is more likely tied to the full climo lifecycle, including some combination
of:

- open
- derived-variable access
- subset/load/detach
- repeated reopen of the same file within one task
- close timing

## Current Code State

Two targeted changes are already in place in `dataset_xr.py`.

1. `_subset_vars_and_load()` now loads synchronously and can return a deep,
   detached in-memory dataset.
2. `_open_climo_dataset()` now uses `xc.open_dataset()` for a concrete
   single-file climatology path and keeps `xc.open_mfdataset()` only for true
   multi-file or globbed climatologies.

These changes were intended to narrow the failure mode, not to claim a fix.

## Relevant Code

- Climo lifecycle wrapper:
  `e3sm_diags/driver/utils/dataset_xr.py`
- Climo open path:
  `e3sm_diags/driver/utils/dataset_xr.py`
- Subset/load/detach logic:
  `e3sm_diags/driver/utils/dataset_xr.py`
- Worker task logging:
  `e3sm_diags/e3sm_diags_driver.py`
- Worker traceback dump setup:
  `e3sm_diags/e3sm_diags_driver.py`
- Focused opener repro:
  `auxiliary_tools/debug/1048-py314-stall-cont/parallel-nersc/min-scripts/open_repro.py`

## Upstream Background And Practical Conclusions

Some useful upstream xarray context matches what this repro is showing.

- In xarray source, the netCDF4 backend still defines
  `NETCDF4_PYTHON_LOCK = combine_locks([NETCDFC_LOCK, HDF5_LOCK])`.
- Source:
  https://github.com/pydata/xarray/blob/main/xarray/backends/netCDF4_.py
- Those base locks are xarray `SerializableLock` objects.
- xarray's own `SerializableLock` docstring says these locks are per-process
  and "will not block concurrent operations between processes."
- Source:
  https://github.com/pydata/xarray/blob/main/xarray/backends/locks.py
- xarray's public `open_dataset()` docs still describe `lock` as the resource
  lock used for safe parallel disk access and say that, by default,
  appropriate locks are chosen based on the active Dask scheduler.
- Docs:
  https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html

That matters here because this repro uses multiprocessing, and the default
netCDF4 backend lock behavior is clearly not sufficient on this Python 3.14
stack.

There is also relevant upstream discussion in xarray PR #10788:

- PR:
  https://github.com/pydata/xarray/pull/10788
- That PR fixed one class of netCDF4 close-related failures.
- But the PR author explicitly called it a partial fix and said that, for the
  `netCDF4` backend, the issue still remains.

## What Still Looks Most Likely

- An xarray/netCDF4 backend lock bug or starvation path under this exact
  multi-process, read-only software stack.
- Repeated reopen churn against the same climatology file within one task.
  The ALBEDO path opens the same staged test climo multiple times for derived
  `ALBEDO`, `LANDFRAC`, `OCNFRAC`, and `landfrac`.
- Interaction between backend file-manager state and the real
  open/load/close lifecycle.
- Shared-filesystem behavior that is still hidden because the current repro
  stages symlinks, not true local copies.

Things that now look less likely:

- plain `fork` inheritance by itself
- single-file misuse of `open_mfdataset()` by itself
- HDF5 file locking alone
- immediate `close()` calls as the main trigger

The last point matters because the main climo path now loads synchronously,
returns a detached dataset, and only then closes the original file-backed
dataset. The observed worker stacks block on lock acquire, not on close.

## Next Steps And Fallbacks

Do not treat `lock=False` as the fix. Use the completed 10-run `lock=False`
result only as evidence that the backend lock path is required for
reproduction.

Best next steps:

1. Build the smallest failing repro that still includes open, derived
   `ALBEDO`, subset/load, and close sequencing.
2. Test true local copies instead of symlinked staging.
3. Reduce repeated opens within one task by opening a climo dataset once per
   filepath, materializing all needed variables, then closing once.

Practical fallback options are narrower than they first appeared.

- Keeping `engine="netcdf4"` on the last known-good Python version remains the
  most straightforward fallback.
