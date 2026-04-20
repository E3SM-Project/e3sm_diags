## Summary

I can still reproduce intermittent stalls on Python 3.14 with `xarray>=2026.01.0`
using the NERSC minimum script (`num_workers=8`). Some runs finish; some
hang.

The current evidence still points to the xarray/netCDF4 lock-acquire path.
The next multiprocessing start-method experiment should be `forkserver`
instead of `fork`.

## Key finding

A targeted `SIGUSR1` on the actual stuck worker removed the attribution
ambiguity.

- Stuck worker: `pid=1313623`
- Last `e3sm_diags` log line from that worker: `Climo backend open_mfdataset start`
- Task context: test `ANN` `ALBEDO` climo file
  `v3.LR.historical_0051_ANN_198501_201412_climo.nc`
- Python stack at the time of the hang was in the xarray/netCDF4 lock path:

```text
xarray/backends/locks.py:66  __enter__
xarray/backends/file_manager.py:181  _optional_lock
xarray/backends/file_manager.py:217  _acquire_with_cache_info
xarray/backends/netCDF4_.py:532  _acquire
xarray/backends/netCDF4_.py:108  get_array
```

This points to the backend open/read path, not a downstream plotting or
regridding phase.

I also re-ran after changing single-file climatologies to use
`xc.open_dataset()` instead of `xc.open_mfdataset()`. The stall still
reproduced.

- Unmatched task in the newer run: `pid=1954444`
- Last completed task from that worker: `HadISST` `DJF`
- Next task that started and never logged `done`: `ceres_ebaf_toa_v4.1` `ANN`
- Last `e3sm_diags` log line from that worker in the newer run:
  `Climo backend open_dataset start`
- Worker traceback sink still showed the same xarray/netCDF4 lock path

So the stall is not limited to single-file use of `open_mfdataset()`.

## Relevant code

- Climo lifecycle wrapper:
  https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/e3sm_diags/driver/utils/dataset_xr.py#L585-L657
- Climo open path:
  https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/e3sm_diags/driver/utils/dataset_xr.py#L664-L760
- Subset/load/detach logic:
  https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/e3sm_diags/driver/utils/dataset_xr.py#L1785-L1835
- Worker task logging:
  https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/e3sm_diags/e3sm_diags_driver.py#L382-L408
- Worker traceback dump setup:
  https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/e3sm_diags/e3sm_diags_driver.py#L427-L477

## Current code state

I made two targeted changes in `dataset_xr.py`.

1. `_subset_vars_and_load()` now loads the subset synchronously and can return
   a deep-copied detached dataset so callers can close the original live
   dataset after load.
2. `_open_climo_dataset()` now detects a concrete single-file climatology path
   and uses `xc.open_dataset()` for that case; it keeps
   `xc.open_mfdataset()` only for true multi-file/globbed climatologies.

The second change is intentional because `open_mfdataset()` still takes the
multi-file machinery path even when the input resolves to one file.

## Focused repro result

I ran a focused concurrent-open repro against the exact `ANN` climo file with
8 worker processes x 20 iterations/worker.

- `xc.open_mfdataset()`: 160/160 successful opens, no timeouts
- `xc.open_dataset()`: 160/160 successful opens, no timeouts
- `xr.open_dataset()`: 160/160 successful opens, no timeouts

Repro script:

- https://github.com/E3SM-Project/e3sm_diags/blob/bugfix/py314-stall-cont/auxiliary_tools/debug/1048-py314-stall-cont/parallel-nersc/min-scripts/open_repro.py#L1-L181

So the simple open-only repro is not sufficient to trigger the hang by itself.
The stall appears to depend on more of the real climo lifecycle around the open
call.

## Current interpretation

- The hang is still strongly correlated with the climo backend open path under
  multiprocessing.
- The targeted stack shows the worker blocked in xarray/netCDF4 lock
  acquisition.
- The newer run shows that this still happens after switching single-file
  climatologies to `xc.open_dataset()`.
- The simplified opener-only repro does not hang, so the trigger is likely the
  full `e3sm_diags` lifecycle around open, derived-variable access,
  subset/load, and close timing.

## Next step

Re-run the ALBEDO repro with `forkserver` instead of `fork`.

- If the stall disappears, inherited post-fork state is likely part of the
  problem.
- If the stall remains, the next focused repro should include derived `ALBEDO`,
  `_subset_vars_and_load(..., detach=True)`, and `ds.close()` sequencing.
