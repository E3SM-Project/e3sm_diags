# Python 3.14 Stall Notes

## Current status

I can still reproduce intermittent stalls on Python 3.14 with
`xarray>=2026.01.0` using the NERSC minimum script with `num_workers=8`.

Current version split:

- `xarray==2025.12.0`: no stall observed
- `xarray>=2026.01.0`: intermittent stall reproduces
- `xarray>=2026.01.0` with `lock=False` in the climo open path: 10/10 runs completed

This currently points to the xarray `netcdf4` backend lock path as the
strongest requirement for reproduction. Because this workflow is read-only and
the climo path eagerly loads and detaches data, `lock=False` is now a realistic
workaround candidate, but it still needs production-scale validation before it
can be treated as safe.

## Key findings

The current evidence points to the xarray/netCDF4 backend lock path as the
strongest current lead. Several obvious reductions still reproduce the stall,
while multiple focused repros do not.

| Test                                                                        | Result           | Interpretation                                                   |
| --------------------------------------------------------------------------- | ---------------- | ---------------------------------------------------------------- |
| `lock=False` in climo open path                                             | no stall ✅      | 10/10 runs completed; default backend lock path appears required |
| plain `xr.open_dataset()` / `xr.open_mfdataset()` instead of xCDAT wrappers | stalls ❌        | not xCDAT-specific                                               |
| single-file `xc.open_dataset()` for concrete climatology files              | stalls ❌        | not just single-file misuse of `open_mfdataset()`                |
| direct `multiprocessing.Pool` instead of `dask.bag`                         | stalls ❌        | `dask.bag` is not the main differentiator                        |
| `forkserver` instead of `fork`                                              | stalls ❌        | not just inherited post-`fork` state                             |
| focused open-only repro                                                     | no stall ✅      | open-only activity is insufficient                               |
| focused ALBEDO lifecycle repro                                              | no stall ✅      | reduced single-file lifecycle is insufficient                    |
| focused ALBEDO lifecycle repro with repeated same-file reopens              | no stall ✅      | repeated same-file reopen alone is insufficient                  |
| true local-copy staging in the minimum script                               | stalls ❌        | not just symlink or staged-filesystem aliasing                   |
| one-open land/sea mask helper in the minimum script                         | stalls ❌        | removed mask alias reopens, but stall still reproduced           |
| `h5netcdf` engine check                                                     | fails to open ❌ | not a viable fallback for this dataset                           |

## What the hang looks like

Some parallel runs finish and some hang.

When the hang occurs, the stuck worker is not spending time in plotting,
regridding, or later output steps. The best worker traceback so far points to
the xarray/netCDF4 backend lock path during climatology open/read work.

A targeted `SIGUSR1` removed ambiguity about which worker was actually stuck.

First clear case:

- stuck worker: `pid=1313623`
- task context: test `ANN` `ALBEDO` climo file
  `v3.LR.historical_0051_ANN_198501_201412_climo.nc`
- last `e3sm_diags` log line from that worker:
  `Climo backend open_mfdataset start`

Representative stack at hang time:

    xarray/backends/locks.py:66  __enter__
    xarray/backends/file_manager.py:181  _optional_lock
    xarray/backends/file_manager.py:217  _acquire_with_cache_info
    xarray/backends/netCDF4_.py:532  _acquire
    xarray/backends/netCDF4_.py:108  get_array

Conclusion: the observed stall is in the climo backend open/read path, inside
the xarray `netcdf4` backend lock/file-manager path.

Latest confirming case after the one-open land/sea mask change:

- stuck worker: `pid=2251775`
- worker traceback file:
  `run-001/prov/worker-stacks/pid-2251775.log`
- representative stack still shows:

  xarray/backends/locks.py:66 **enter**
  xarray/backends/locks.py:229 **enter**
  xarray/backends/file*manager.py:181 \_optional_lock
  xarray/backends/file_manager.py:217 \_acquire_with_cache_info
  xarray/backends/netCDF4*.py:532 _acquire
  xarray/backends/netCDF4_.py:108 get_array

Conclusion: even after removing the repeated mask alias reopen sequence, the
actual stuck worker is still blocking in the same xarray/netCDF4 backend lock
path.

## Experiment details

### Diagnostic evidence

#### 1. `lock=False` in climo open path

I re-ran the same looped ALBEDO repro with `lock=False` passed through the
xarray/xCDAT open calls.

Result:

- all 10 `lock=False` runs completed without stalling
- earlier default-lock repros usually stalled by run 2

Conclusion: the default xarray `netcdf4` backend lock path appears to be
required for reproduction in this workflow. `lock=False` should be treated as
diagnostic evidence only.

### Still reproduces with reductions

#### 2. Plain `xr.open_dataset()` / `xr.open_mfdataset()` instead of xCDAT wrappers

I replaced `xc.open_dataset()` and `xc.open_mfdataset()` with plain
`xr.open_dataset()` and `xr.open_mfdataset()` in the climatology open path.

Result:

- the stall still reproduced
- the worker traceback still pointed to the same xarray/netCDF4 lock-acquire
  path during climatology open/read work

Conclusion: xCDAT does not appear to be the essential differentiator for this
stall.

#### 3. Single-file `xc.open_dataset()` for concrete climatology files

I changed single concrete climatology paths to use `xc.open_dataset()` and kept
`xc.open_mfdataset()` only for true multi-file or globbed inputs.

Result:

- the stall still reproduced
- the later stuck worker showed the last log line
  `Climo backend open_dataset start`
- the worker traceback still showed the same xarray/netCDF4 lock path

Conclusion: this is not just a single-file misuse of `open_mfdataset()`.

#### 4. `forkserver` instead of `fork`

I re-ran with `forkserver`.

Result:

- the stall still reproduced
- the stuck worker traceback sink was `pid=2254799`
- that traceback still showed the same xarray/netCDF4 lock-acquire path

Conclusion: plain post-`fork` inherited state is not sufficient to explain the
stall, although a broader process-model interaction is still possible.

#### 4b. Direct `multiprocessing.Pool` instead of `dask.bag`

I replaced the `dask.bag` process scheduler path in `e3sm_diags_driver.py`
with a direct `multiprocessing.Pool(...).map(...)` path while keeping the same
task list, worker count, fork start method, and worker traceback sink.

Result:

- the run still stalled
- the stuck worker traceback still showed the same xarray/netCDF4 backend lock
  path
- a worker-side `SIGUSR1` dump for `pid=1385846` again showed:

  xarray/backends/locks.py:66 **enter**
  xarray/backends/locks.py:229 **enter**
  xarray/backends/file*manager.py:181 \_optional_lock
  xarray/backends/file_manager.py:217 \_acquire_with_cache_info
  xarray/backends/netCDF4*.py:532 _acquire
  xarray/backends/netCDF4_.py:108 get_array

Conclusion: `dask.bag` is not the main differentiator for this stall.
Replacing it with direct standard-library multiprocessing simplifies the
execution path, but it does not change the core failure signature.

### Did not reproduce in focused repros

#### 5. Focused open-only repro

I ran a focused concurrent-open repro against the exact `ANN` climo file with
8 worker processes and 20 iterations per worker.

Result:

- `xc.open_mfdataset()`: 160/160 successful opens, no timeouts
- `xc.open_dataset()`: 160/160 successful opens, no timeouts
- `xr.open_dataset()`: 160/160 successful opens, no timeouts

Conclusion: simple concurrent open-only activity is not enough to trigger the
hang.

#### 6. Focused ALBEDO lifecycle repro

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

#### 7. Focused ALBEDO lifecycle repro with repeated same-file reopens

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

#### 8. `h5netcdf` engine check

I tested whether `h5netcdf` could serve as a backend workaround for the same
staged climatology file.

Result:

- the open failed immediately with
  `OSError: Unable to synchronously open file (file signature not found)`
- `ncdump -k` on the same file reported `64-bit offset`

Conclusion: these staged climo files are not HDF5-backed netCDF4 files, so
`h5netcdf` is not currently a usable fallback for this dataset.

#### 9. True local-copy staging in the minimum script

I re-ran the NERSC minimum script using true local copies instead of symlinked
staging for the model climatology inputs.

Result:

- the run still stalled
- the logs showed the same repeated same-file climo open pattern for
  `ALBEDO`, `LANDFRAC`, `OCNFRAC`, and `landfrac`
- the traceback captured from the notebook process only showed the main thread
  blocked waiting for Dask completion, not a contradictory worker-side cause

Conclusion: the stall is not explained by symlink staging or a simple staged
filesystem alias effect. This shifts weight away from the staging mechanism and
toward workflow shape, especially repeated opens across Dask tasks.

#### 10. One-open land/sea mask helper in the minimum script

I changed the land/sea mask helper so one task opens the test climo file once
for mask extraction instead of reopening it separately for `LANDFRAC`,
`OCNFRAC`, and `landfrac` alias checks.

Result:

- the repeated same-file mask alias reopen sequence disappeared from the logs
- the test climo file was still opened a second time within the task after the
  main `ALBEDO` load/derive path
- the run still stalled after the `DJF` task completed
- the notebook `SIGUSR1` traceback again only showed the parent process blocked
  in Dask wait, not the stuck worker-side frame
- a worker-side `SIGUSR1` dump for `pid=2251775` confirmed the actual stuck
  worker was still blocked in `xarray.backends.locks` /
  `xarray.backends.file_manager` /
  `xarray.backends.netCDF4_._acquire`

Conclusion: the repeated mask alias reopens were real churn, but eliminating
them was not sufficient to remove the stall. The remaining in-task duplicate
open on the same climo file is now the more relevant low-level target if the
current parallel structure must be preserved.

## Current interpretation

### Most likely

- an xarray `netcdf4` backend lock bug, starvation path, or lock-related
  regression under this exact multi-process read-only software stack
- interaction between backend file-manager state and a broader workflow than
  the focused staged-test-file repros currently exercise
- workflow-specific interaction involving additional files, reference-data
  access, repeated same-file opens across parallel tasks, or later driver steps
- repeated open/read of the same climatology file inside one task, even after
  the mask alias reopen sequence was reduced to a single additional open

### Less likely

- xCDAT-specific wrapper behavior by itself
- `dask.bag` by itself
- single-file misuse of `open_mfdataset()` by itself
- plain `fork` inheritance by itself
- immediate `close()` calls as the main trigger
- simple concurrent open activity by itself
- symlink staging or a simple local-copy versus symlink difference by itself
- the land/sea mask alias reopen sequence by itself

That last point matters because the main climo path now loads synchronously,
returns a detached dataset, and only then closes the original file-backed
dataset. The observed worker stacks block on lock acquire, not on close.

## Practical implication

The trigger is probably not just parallel open activity by itself.

The remaining gap is between the focused repros that succeed and the full
driver workflow that still stalls. The missing trigger is more likely to be in
workflow composition, such as additional file opens, reference-data access,
file-manager interactions across a longer task lifecycle, or the remaining
duplicate same-file open inside a task after the mask helper reduction.

## Current code state

Two targeted changes are already in place in `dataset_xr.py`:

1. `_subset_vars_and_load()` now loads synchronously and can return a deep,
   detached in-memory dataset.
2. `_get_default_land_sea_mask_dataset()` now loads, deep-copies, and closes
   the fallback mask dataset before returning it.
3. `_get_land_sea_mask_dataset()` now opens the seasonal climo file once for
   mask extraction instead of reopening it separately for mask alias checks.

These changes were made to narrow the failure mode, not to claim a fix.

The multiprocessing driver path in `e3sm_diags_driver.py` now uses direct
`multiprocessing.Pool` instead of `dask.bag`, which simplifies orchestration
but did not change the stall behavior.

## Relevant code

- `e3sm_diags/driver/utils/dataset_xr.py`
  - climo lifecycle wrapper
  - climo open path
  - subset/load/detach logic
- `e3sm_diags/e3sm_diags_driver.py`
  - direct multiprocessing task execution
  - worker task logging
  - worker traceback dump setup
- `auxiliary_tools/debug/1048-py314-stall-cont/parallel-nersc/min-scripts/open_repro.py`
  - focused opener repro

## Upstream background

Relevant xarray background lines up with the current repro.

The `netCDF4` backend in xarray uses a combined backend lock:

- `NETCDF4_PYTHON_LOCK = combine_locks([NETCDFC_LOCK, HDF5_LOCK])`
- Source: <https://github.com/pydata/xarray/blob/main/xarray/backends/netCDF4_.py>

Those underlying locks are xarray `SerializableLock` objects. In xarray source,
the `SerializableLock` docstring says they are per-process and "will not block
concurrent operations between processes".

- Source: <https://github.com/pydata/xarray/blob/main/xarray/backends/locks.py>

xarray's public `open_dataset()` documentation also describes `lock` as the
resource lock used for safe parallel disk access, with defaults chosen based on
the active Dask scheduler.

- Docs: <https://docs.xarray.dev/en/stable/generated/xarray.open_dataset.html>

This matters because the repro uses multiprocessing, and the observed behavior
suggests the default `netCDF4` backend locking path is implicated on this
Python 3.14 stack.

There is also relevant upstream discussion in xarray PR `#10788`:

- PR: <https://github.com/pydata/xarray/pull/10788>
- it fixed one class of `netCDF4` close-related failures
- the PR author explicitly described it as a partial fix
- for the `netCDF4` backend, the broader issue still remained

There is also older but directly relevant xarray discussion in issue `#824`
("Disable lock=True in open_mfdataset when reading netCDF3 files"):

- `shoyer` explicitly said the no-lock concern applies to `netCDF4`/HDF5, and
  that with `netCDF3` you do not need the lock when reading data
- the same comment says `lock=False` should not be used blindly for all cases;
  instead `lock=None` should choose something smart based on file type/backend

This matters here because the stalled climatology file is `netCDF3`, the
workload is read-only, and the worker is still blocking in the xarray
`netCDF4_` backend lock path. That makes `lock=False` a more principled
workaround candidate for this specific path, and it also strengthens the case
that the current behavior may be a regression or an uncovered backend path.

## Not yet shown

- not yet shown that this is an xarray-only issue independent of E3SM-Diags
  workflow composition
- not yet shown that `lock=False` is numerically safe for production use at
  full scale

## Next steps

1. run the full-scale production case with `lock=False` at least 3 times
2. compare metrics, JSON outputs, and any saved NetCDF products against a
   trusted baseline (`xarray==2025.12.0` if available, otherwise a known-good
   serial or non-stalling run)
3. compare the repeated `lock=False` runs against each other to confirm they
   are numerically stable, not just stall-free
4. if those checks pass, treat `lock=False` as a temporary workaround candidate
   for `xarray>=2026.01.0` while preparing an upstream report
5. add a minimal version-boundary summary to any upstream report:
   `2025.12.0` works, `2026.01.0+` stalls, `lock=False` removes the stall
6. in any upstream report, include that xarray issue `#824` says `netCDF3`
   reads should not need the lock, while this read-only `netCDF3` path still
   appears to block in lock acquisition
