# Python 3.14 stall investigation

## Current status

I can still reproduce intermittent stalls on Python 3.14 with `xarray>=2026.01.0` using the NERSC minimum script with `num_workers=8`.

Version boundary:

- `xarray==2025.12.0`: no stall observed
- `xarray>=2026.01.0`: intermittent stall reproduces
- `xarray>=2026.01.0` with `lock=False` in the climo open path: 10/10 runs completed

Current best lead: the xarray `netcdf4` backend lock path during climo open/read work.

`lock=False` is now a plausible workaround candidate for this read-only `NetCDF3` path, but it still needs production-scale validation.

## Why xarray is the leading suspect

- The version split is clean: `2025.12.0` works, `2026.01.0+` stalls.
- Worker traces consistently point to the same xarray/netCDF4 lock-acquire path.
- The stall persists after removing xCDAT wrappers and after a temporary direct
  `multiprocessing.Pool` replacement for `dask.bag`.
- `lock=False` removes the stall in the current repro.

## Key findings

| Test                                                                        | Result           | Takeaway                                                         |
| --------------------------------------------------------------------------- | ---------------- | ---------------------------------------------------------------- |
| `lock=False` in climo open path                                             | no stall ✅      | 10/10 runs completed; default backend lock path appears required |
| plain `xr.open_dataset()` / `xr.open_mfdataset()` instead of xCDAT wrappers | stalls ❌        | not xCDAT-specific                                               |
| single-file `xc.open_dataset()` for concrete climatology files              | stalls ❌        | not just single-file misuse of `open_mfdataset()`                |
| temporary direct `multiprocessing.Pool` instead of `dask.bag`               | stalls ❌        | `dask.bag` is not the main differentiator                        |
| `forkserver` instead of `fork`                                              | stalls ❌        | not just inherited post-fork state                               |
| focused open-only repro                                                     | no stall ✅      | open-only activity is insufficient                               |
| focused ALBEDO lifecycle repro                                              | no stall ✅      | reduced single-file lifecycle is insufficient                    |
| focused ALBEDO lifecycle repro with repeated same-file reopens              | no stall ✅      | same-file reopen alone is insufficient                           |
| true local-copy staging in the minimum script                               | stalls ❌        | not explained by symlink or staged-filesystem aliasing           |
| one-open land/sea mask helper in the minimum script                         | stalls ❌        | removing mask alias reopens was not enough                       |
| `h5netcdf` engine check                                                     | fails to open ❌ | not a viable fallback for this dataset                           |

## What the hang looks like

Some parallel runs finish and some hang.

When the hang occurs, the stuck worker is not in plotting, regridding, or later output steps. The best worker traceback so far points to the xarray/netCDF4 backend lock path during climo open/read work.

A targeted `SIGUSR1` dump removed ambiguity about which worker was actually stuck.

### First clear stuck case

- stuck worker: `pid=1313623`
- task context: `ANN` `ALBEDO`
- file: `v3.LR.historical_0051_ANN_198501_201412_climo.nc`
- file format: `NetCDF3`
- last worker log line: `Climo backend open_mfdataset start`

Representative stack:

    xarray/backends/locks.py:66              __enter__
    xarray/backends/file_manager.py:181      _optional_lock
    xarray/backends/file_manager.py:217      _acquire_with_cache_info
    xarray/backends/netCDF4_.py:532          _acquire
    xarray/backends/netCDF4_.py:108          get_array

### Confirming case after mask-helper reduction

- stuck worker: `pid=2251775`
- traceback file: `run-001/prov/worker-stacks/pid-2251775.log`

Representative stack:

    xarray/backends/locks.py:66              __enter__
    xarray/backends/locks.py:229             __enter__
    xarray/backends/file_manager.py:181      _optional_lock
    xarray/backends/file_manager.py:217      _acquire_with_cache_info
    xarray/backends/netCDF4_.py:532          _acquire
    xarray/backends/netCDF4_.py:108          get_array

Conclusion: even after reducing repeated mask reopens, the actual stuck worker is still blocking in the same xarray/netCDF4 backend lock path.

## Experiment summary

### Reproduces with default locking, disappears with `lock=False`

#### `lock=False` in climo open path

I re-ran the looped ALBEDO repro with `lock=False` passed through the xarray/xCDAT open calls.

Result:

- all 10 `lock=False` runs completed
- earlier default-lock repros usually stalled by run 2

Conclusion: the default xarray `netcdf4` backend lock path appears required for reproduction in this workflow.

### Still reproduces after reductions

#### Plain xarray open calls instead of xCDAT wrappers

Replaced `xc.open_dataset()` and `xc.open_mfdataset()` with plain xarray calls.

Result:

- stall still reproduced
- stuck worker traceback still pointed to the same xarray/netCDF4 lock-acquire path

Conclusion: xCDAT is not the essential differentiator.

#### Single-file `xc.open_dataset()` for concrete climatology files

Used `xc.open_dataset()` for concrete single-file inputs and reserved `xc.open_mfdataset()` for true multi-file or globbed inputs.

Result:

- stall still reproduced
- later stuck worker showed `Climo backend open_dataset start`
- traceback still showed the same xarray/netCDF4 lock path

Conclusion: not just single-file misuse of `open_mfdataset()`.

#### `forkserver` instead of `fork`

Re-ran with `forkserver`.

Result:

- stall still reproduced
- traceback still showed the same xarray/netCDF4 lock-acquire path

Conclusion: plain post-fork inherited state is not enough to explain the stall.

#### Direct `multiprocessing.Pool` instead of `dask.bag`

Replaced the `dask.bag` process scheduler path with direct `multiprocessing.Pool(...).map(...)`.

Result:

- run still stalled
- stuck worker traceback still showed the same xarray/netCDF4 backend lock path

Conclusion: `dask.bag` is not the main differentiator.

### Focused repros that did not reproduce

#### Focused open-only repro

Ran a concurrent-open repro against the exact `ANN` climo file with 8 workers and 20 iterations per worker.

Result:

- `xc.open_mfdataset()`: 160/160 successful opens
- `xc.open_dataset()`: 160/160 successful opens
- `xr.open_dataset()`: 160/160 successful opens

Conclusion: simple concurrent open-only activity is insufficient.

#### Focused ALBEDO lifecycle repro

Ran a reduced repro closer to the real ALBEDO path: open, derive `ALBEDO`, subset/load/detach, and close.

Result:

- 8 workers x 20 iterations completed
- 160/160 iterations succeeded
- no timeouts or worker errors
- real source tuple used: `("SOLIN", "FSNTOA")`

Conclusion: the single-file ALBEDO derive/load/detach/close lifecycle is insufficient.

#### Focused ALBEDO lifecycle repro with repeated same-file reopens

Extended the focused repro so each iteration also reopened the same staged climo file for fraction variables after the ALBEDO lifecycle.

Result:

- 8 workers x 20 iterations completed
- 160/160 iterations succeeded
- no timeouts or worker errors

Conclusion: same-file reopen alone is insufficient.

#### `h5netcdf` engine check

Tested whether `h5netcdf` could be used as a backend workaround.

Result:

- open failed immediately with `OSError: Unable to synchronously open file (file signature not found)`
- `ncdump -k` reported `64-bit offset`

Conclusion: `h5netcdf` is not a viable fallback for this dataset.

### Full minimum script still reproduces

#### True local-copy staging

Re-ran the NERSC minimum script using true local copies instead of symlinked staging.

Result:

- run still stalled
- logs still showed repeated same-file climo opens for `ALBEDO`, `LANDFRAC`, `OCNFRAC`, and `landfrac`

Conclusion: not explained by symlink staging or filesystem aliasing.

#### One-open land/sea mask helper

Changed the land/sea mask helper so one task opens the test climo file once for mask extraction instead of reopening it separately for `LANDFRAC`, `OCNFRAC`, and `landfrac`.

Result:

- repeated mask alias reopen sequence disappeared from the logs
- the file was still opened a second time later within the task
- run still stalled
- worker-side `SIGUSR1` again confirmed blocking in `xarray.backends.locks` / `file_manager` / `netCDF4_._acquire`

Conclusion: removing the mask alias reopen sequence was not enough. The remaining duplicate same-file open inside a task is now the more relevant low-level target.

## Current interpretation

Most likely:

- a lock-related regression or starvation path in the xarray `netcdf4` backend under this Python 3.14 stack
- interaction between backend file-manager state and the broader full-driver workflow
- a workflow-composition trigger not present in the focused repros, such as additional file opens, reference-data access, longer task lifecycles, or repeated same-file open/read within a task

Less likely:

- xCDAT-specific behavior alone
- `dask.bag` alone
- single-file misuse of `open_mfdataset()` alone
- plain `fork` inheritance alone
- simple concurrent open activity alone
- symlink staging or aliasing alone
- the mask alias reopen sequence alone
- immediate `close()` as the main trigger

## Current code state

These targeted changes are already in place in `dataset_xr.py`:

1. `_subset_vars_and_load()` now loads synchronously and can return a deep, detached in-memory dataset.
2. `_get_default_land_sea_mask_dataset()` now loads, deep-copies, and closes the fallback mask dataset before returning it.
3. `_get_land_sea_mask_dataset()` now opens the seasonal climo file once for mask extraction instead of reopening it separately for mask alias checks.

These changes narrow the failure mode. They do not fix it.

In `e3sm_diags_driver.py`, the code remains on the original `dask.bag`
multiprocessing path to keep this PR scoped. A temporary direct
`multiprocessing.Pool(...).map(...)` experiment also reproduced the stall and
was reverted.

## Upstream background

Relevant xarray details line up with the current repro:

- the xarray `netCDF4` backend uses a combined backend lock:
  `NETCDF4_PYTHON_LOCK = combine_locks([NETCDFC_LOCK, HDF5_LOCK])`
- those underlying locks are xarray `SerializableLock` objects
- the `SerializableLock` docstring says they are per-process and "will not block concurrent operations between processes"
- xarray `open_dataset()` docs describe `lock` as the resource lock for safe parallel disk access, with defaults chosen based on scheduler

Why that matters here:

- the repro uses multiprocessing
- the stuck file is `NetCDF3`
- the workload is read-only
- worker traces still block in the xarray `netCDF4_` backend lock path

Related upstream discussion:

- xarray PR `#10788`: fixed one class of `netCDF4` close-related failures, but was described as only a partial fix
- xarray issue `#824`: discusses that netCDF3 reads should not need the lock, while also noting `lock=False` should not be applied blindly in all cases

This strengthens the case that `lock=False` is a principled workaround candidate for this specific read-only `NetCDF3` path, and that the current behavior may reflect a regression or uncovered backend path.

## Workaround scope and alternatives

`lock=False` is acceptable only for verified, read-only, physically NetCDF3 climo files. The check should require every file to report a NetCDF3 data model and must exclude `NETCDF4_CLASSIC`, which is still HDF5-backed.

Keep xarray's default lock behavior for NetCDF4, `NETCDF4_CLASSIC`, mixed, unreadable, or unknown inputs.

Alternatives:

- **Default lock everywhere:** safest general behavior, but preserves the stall.
- **`engine="scipy"` for NetCDF3:** cleaner backend boundary, but a larger behavior change.
- **Disable lock globally:** rejected; unsafe for NetCDF4/HDF5 parallel reads.
- **Eager load:** may avoid lazy lock paths, but increases memory and does not target the backend issue.

## Not yet shown

- that this is an xarray-only issue independent of E3SM-Diags workflow composition
- that `lock=False` is numerically safe for full production use

## Next steps

1. Run the full production case with `lock=False` at least 3 times.
2. Compare metrics, JSON outputs, and saved NetCDF products against a trusted baseline, ideally `xarray==2025.12.0`.
3. Compare the repeated `lock=False` runs against each other for numerical stability.
4. If those checks pass, treat `lock=False` as a temporary workaround candidate for `xarray>=2026.01.0`.
5. Prepare an upstream report with the minimal version boundary:
   - `2025.12.0` works
   - `2026.01.0+` stalls
   - `lock=False` removes the stall
6. In that report, include that xarray issue `#824` indicates netCDF3 reads should not need the lock, while this read-only netCDF3 path still appears to block in lock acquisition.