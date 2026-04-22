# Python 3.14 Stall Notes

## Current status

I can still reproduce intermittent stalls on Python 3.14 with
`xarray>=2026.01.0` using the NERSC minimum script with `num_workers=8`.

The strongest current evidence points to the xarray/netCDF4 backend
lock-acquire path during climatology open/read work.

## Key findings

### Strongest positive signal

A 10-run diagnostic using `lock=False` in the climo open path completed without
stalling, while earlier `lock=True` repros typically stalled by run 2.

This is strong evidence that the xarray/netCDF4 backend lock path is required
for reproduction.

Do not treat `lock=False` as the fix. Treat it only as diagnostic evidence.

### Strongest negative signals

The stall still reproduces after each of the following reductions:

- switching single-file climatologies from `xc.open_mfdataset()` to
  `xc.open_dataset()`
- switching xCDAT open calls to plain `xr.open_dataset()` /
  `xr.open_mfdataset()`
- switching multiprocessing from `fork` to `forkserver`

This makes the following explanations less likely as primary causes:

- single-file misuse of `open_mfdataset()`
- xCDAT wrapper behavior by itself
- plain post-`fork` inherited state by itself

### Result matrix

| Test                                                                        | Result                  | Interpretation                                    |
| --------------------------------------------------------------------------- | ----------------------- | ------------------------------------------------- |
| single-file `xc.open_dataset()` for concrete climatology files              | still stalls            | not just single-file misuse of `open_mfdataset()` |
| plain `xr.open_dataset()` / `xr.open_mfdataset()` instead of xCDAT wrappers | still stalls            | not xCDAT-specific                                |
| `forkserver` instead of `fork`                                              | still stalls            | not just post-`fork` inherited state              |
| `lock=False` in climo open path                                             | 10/10 runs completed    | backend lock path likely required                 |
| focused open-only repro                                                     | no stall                | open-only activity is insufficient                |
| focused ALBEDO lifecycle repro                                              | no stall                | reduced single-file lifecycle is insufficient     |
| focused ALBEDO lifecycle repro with repeated same-file reopens              | no stall                | repeated same-file reopen alone is insufficient   |
| `h5netcdf` engine check                                                     | unusable on staged file | not a viable fallback for this dataset            |

## What the hang looks like

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

## Reproduction details

### Reproduces

#### 1. Targeted stuck-worker traceback

A targeted `SIGUSR1` removed ambiguity about which worker was actually stuck.

First clear case:

- stuck worker: `pid=1313623`
- task context: test `ANN` `ALBEDO` climo file
  `v3.LR.historical_0051_ANN_198501_201412_climo.nc`
- last `e3sm_diags` log line from that worker:
  `Climo backend open_mfdataset start`

Conclusion: the hang is tied to the climo backend open/read path.

#### 2. Single-file `open_dataset()` instead of `open_mfdataset()`

I changed single concrete climatology paths to use `xc.open_dataset()` and kept
`xc.open_mfdataset()` only for true multi-file or globbed inputs.

Result:

- the stall still reproduced
- the later stuck worker showed the last log line
  `Climo backend open_dataset start`
- the worker traceback still showed the same xarray/netCDF4 lock path

Conclusion: this is not just a single-file misuse of `open_mfdataset()`.

#### 3. `forkserver` instead of `fork`

I re-ran with `forkserver`.

Result:

- the stall still reproduced
- the stuck worker traceback sink was `pid=2254799`
- that traceback still showed the same xarray/netCDF4 lock-acquire path

Conclusion: plain post-`fork` inherited state is not required for reproduction.

#### 4. Plain xarray open calls instead of xCDAT open calls

I replaced `xc.open_dataset()` and `xc.open_mfdataset()` with plain
`xr.open_dataset()` and `xr.open_mfdataset()` in the climatology open path.

Result:

- the stall still reproduced
- the worker traceback still pointed to the same xarray/netCDF4 lock-acquire
  path during climatology open/read work

Conclusion: xCDAT does not appear to be the essential differentiator for this
stall.

#### 5. `lock=False` as a diagnostic

I re-ran the same looped ALBEDO repro with `lock=False` passed through the
xarray/xCDAT open calls.

Result:

- all 10 `lock=False` runs completed without stalling
- earlier `lock=True` repros usually stalled by run 2

Conclusion: the xarray/netCDF4 backend lock path appears to be required for
reproduction.

### Did not reproduce in focused repros

#### 6. Focused open-only repro

I ran a focused concurrent-open repro against the exact `ANN` climo file with
8 worker processes and 20 iterations per worker.

Result:

- `xc.open_mfdataset()`: 160/160 successful opens, no timeouts
- `xc.open_dataset()`: 160/160 successful opens, no timeouts
- `xr.open_dataset()`: 160/160 successful opens, no timeouts

Conclusion: simple concurrent open-only activity is not enough to trigger the
hang.

#### 7. Focused ALBEDO lifecycle repro

I ran a middle-sized repro that exercises the real ALBEDO path more closely:
open the climo file, derive `ALBEDO` from the actual source variables,
subset/load/detach, and close, under multiprocessing.

Result:

- 8 workers x 20 iterations per worker completed successfully
- 160/160 iterations completed
- no worker timeouts
- no worker errors
- the derived ALBEDO path used the real source tuple `("SOLIN", "FSNTOA")`

Conclusion: the single-file ALBEDO derive/load/detach/close lifecycle by
itself is not sufficient to trigger the stall.

#### 8. Focused ALBEDO lifecycle repro with repeated same-file reopens

I then extended the focused ALBEDO repro so each iteration also reopened the
same staged climo file for available fraction variables after the ALBEDO
lifecycle.

Result:

- 8 workers x 20 iterations per worker completed successfully
- 160/160 iterations completed
- no worker timeouts
- no worker errors
- the repeated reopen path materialized both `LANDFRAC` and `OCNFRAC`

Conclusion: repeated same-file reopens for the staged test climo alone are not
sufficient to trigger the stall.

#### 9. `h5netcdf` engine check

I tested whether `h5netcdf` could serve as a backend workaround for the same
staged climatology file.

Result:

- the open failed immediately with
  `OSError: Unable to synchronously open file (file signature not found)`
- `ncdump -k` on the same file reported `64-bit offset`

Conclusion: these staged climo files are not HDF5-backed netCDF4 files, so
`h5netcdf` is not currently a usable fallback for this dataset.

## Current interpretation

### Most likely

- an xarray/netCDF4 backend lock bug or starvation path under this exact
  multi-process, read-only software stack
- interaction between backend file-manager state and a broader workflow than
  the focused staged-test-file repros currently exercise
- workflow-specific interaction involving additional files, reference-data
  access, or later driver steps
- shared-filesystem behavior that is still hidden because the current repro
  stages symlinks, not true local copies

### Less likely

- plain `fork` inheritance by itself
- single-file misuse of `open_mfdataset()` by itself
- xCDAT-specific wrapper behavior by itself
- HDF5 file locking alone
- immediate `close()` calls as the main trigger

That last point matters because the main climo path now loads synchronously,
returns a detached dataset, and only then closes the original file-backed
dataset. The observed worker stacks block on lock acquire, not on close.

## Practical implication

The trigger is probably not just parallel open activity by itself.

The remaining gap is between the focused repros that succeed and the full
driver workflow that still stalls. The missing trigger is more likely to be in
workflow composition, such as additional file opens, reference-data access,
file-manager interactions across a longer task lifecycle, or shared-filesystem
effects that the staged symlink repro still hides.

## Current code state

Two targeted changes are already in place in `dataset_xr.py`:

1. `_subset_vars_and_load()` now loads synchronously and can return a deep,
   detached in-memory dataset.
2. `_open_climo_dataset()` now uses `xc.open_dataset()` for a concrete
   single-file climatology path and keeps `xc.open_mfdataset()` only for true
   multi-file or globbed climatologies.

These changes were made to narrow the failure mode, not to claim a fix.

## Relevant code

- `e3sm_diags/driver/utils/dataset_xr.py`
  - climo lifecycle wrapper
  - climo open path
  - subset/load/detach logic
- `e3sm_diags/e3sm_diags_driver.py`
  - worker task logging
  - worker traceback dump setup
- `auxiliary_tools/debug/1048-py314-stall-cont/parallel-nersc/min-scripts/open_repro.py`
  - focused opener repro

## Upstream background

Relevant xarray background is consistent with what this repro shows:

- in xarray source, the netCDF4 backend defines
  `NETCDF4_PYTHON_LOCK = combine_locks([NETCDFC_LOCK, HDF5_LOCK])`
  Source:
  <https://github.com/pydata/xarray/blob/main/xarray/backends/netCDF4_.py>
- those base locks are xarray `SerializableLock` objects
- xarray's `SerializableLock` docstring says these locks are per-process and
  "will not block concurrent operations between processes"
  Source:
  <https://github.com/pydata/xarray/blob/main/xarray/backends/locks.py>
- xarray's public `open_dataset()` docs describe `lock` as the resource lock
  used for safe parallel disk access, with defaults chosen based on the active
  Dask scheduler
  Docs:
  <https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html>

This matters because the repro uses multiprocessing, and the default netCDF4
backend lock behavior is not sufficient on this Python 3.14 stack.

There is also relevant upstream discussion in xarray PR #10788:

- PR:
  <https://github.com/pydata/xarray/pull/10788>
- that PR fixed one class of netCDF4 close-related failures
- the PR author explicitly described it as a partial fix
- for the `netCDF4` backend, the broader issue still remained

## Next steps

1. test true local copies instead of symlinked staging
2. extend the focused repro toward the next missing pieces in the full driver
   workflow, especially reference-data access or any additional file opens in
   `lat_lon_driver`
3. reduce repeated opens within one task by opening a climo dataset once per
   filepath, materializing all needed variables, then closing once
