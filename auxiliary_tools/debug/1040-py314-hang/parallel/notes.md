# Python 3.14 Parallel Hang Notes

## Summary

- The original long `e3sm_diags` run was stable in serial, hung early with `forkserver`, and stalled later with `fork`.
- A short repro was built in [qa.py](/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/1040-py314-hang/parallel/qa.py) so the problem could be reproduced quickly with the `polar` set and variable `T`.
- The short repro showed that `Z3` is not the unique trigger. The problem also appears on `T`.
- The problem appears under both `forkserver` and `spawn`, so it is not limited to one multiprocessing start method.

## Key Findings

- The first clearly isolated stall happened while loading the land/sea mask, specifically on `test_ds.get_climo_dataset("LANDFRAC", "ANN")`.
- After bypassing the land/sea mask path, a similar stall appeared while opening the main climatology file for `test_ds.get_climo_dataset("T", "ANN")`.
- The hang was narrowed to the shared climatology open path, not to variable derivation, regridding, or one special variable.
- The stall occurred before each of these open calls could return:
  - `xc.open_mfdataset(...)`
  - `xc.open_dataset(...)`
  - `xr.open_dataset(filepath, ...)`

## Strongest Result

- The raw netCDF4-backed open path does not reproduce the stall.
- With `RAW_NETCDF4_CLIMO_OPEN=1`, the narrowed repro completed successfully with `num_workers=1`.
- With `RAW_NETCDF4_CLIMO_OPEN=1`, the narrowed repro also completed successfully with `num_workers=8`.
- This is the strongest evidence so far that the issue is in the normal filepath-based xCDAT/xarray open path, or in related backend/file-manager behavior, rather than in the NetCDF file contents alone.

## Current Hypothesis

- Repeated climatology file opens are eventually wedging in the standard xCDAT/xarray open path.
- A resource-lifecycle problem is plausible because there are no explicit `.close()` calls in the traced `dataset_xr.py` / `io.py` / `polar_driver.py` path.
- Earlier runs also showed relatively high open file descriptor counts, which is consistent with the possibility of retained file-manager or backend state.

## Current Workaround

- Use `RAW_NETCDF4_CLIMO_OPEN=1`.
- This is the best current workaround because it allowed the narrowed parallel repro to finish successfully at both low and higher worker counts.

## Next Step

- Keep `RAW_NETCDF4_CLIMO_OPEN=1` as the working path for now.
- If we want to pursue a deeper fix later, the next best follow-up is to test whether explicit close and cleanup on the standard open path can make the normal xCDAT/xarray route reliable again.
