## Forkserver Debug + Artifact Staging Plan

### Summary
Focus the investigation on the `forkserver` route only and treat `fork` as comparison data, not a supported target for now. Before any code changes, stage the current investigation artifacts into `auxiliary_tools/debug/1040-py314-hang/parallel/artifacts/` so the bash launcher and both failing logs are preserved alongside the local debug driver.

### Key Changes
- Copy these three external artifacts into `auxiliary_tools/debug/1040-py314-hang/parallel/artifacts/` with original basenames preserved:
  - `/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.py314_tom_branch.bash`
  - `/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.tom_branch_multi.o1193363`
  - `/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.tom_branch_multi_fork.o1193430`
- Update the local forkserver debug path so the active repro uses one explicit multiprocessing context and logs it at startup. Prefer a single source of truth in the local debug driver or launcher, not mixed script defaults.
- Add forkserver-only hang instrumentation around the exact pre-stall region:
  - `_get_xarray_datasets()` entry/exit and land-sea-mask fetch in `e3sm_diags/driver/utils/io.py`
  - `_get_land_sea_mask()`, `_get_land_sea_mask_dataset()`, `_get_default_land_sea_mask_dataset()`, `_open_climo_dataset()`, `_subset_vars_and_load()`, and `_get_global_attr_from_climo_dataset()` in `e3sm_diags/driver/utils/dataset_xr.py`
  - `_set_name_yrs_attrs()` call sites already present in `core_parameter.py`; keep those and add worker/process identifiers to all hang logs
  - `polar_driver.run_diag()` around `_get_xarray_datasets()`, `_set_name_yrs_attrs()`, and before `Selected pressure level(s)`
- Include per-log-line process metadata in debug output: PID, parent PID, variable, season, set name, and elapsed time. Add FD count where already supported in `dataset_xr.py`.
- Add a small forkserver runtime guard in the debug path to print the effective worker count and multiprocessing context at process start so future runs cannot be misread from stale scripts.
- Run the forkserver repro in a narrow matrix to isolate the failure trigger:
  - `num_workers=1` with forkserver
  - `num_workers=2` with forkserver
  - `num_workers=4` with forkserver
  - If stable, `num_workers=8` with forkserver
- Use a reduced diagnostic scope first, targeting the observed pre-stall path:
  - Start with `polar`
  - Limit variables to `T` and `Z3`
  - Keep the same seasons/regions that appear near the failures
- If the stall still occurs, use the new logs to classify it into one of four buckets:
  - Child dies during dataset open/load
  - Child hangs during land-sea-mask fetch
  - Child hangs during `get_name_yrs_attr()`
  - Child hangs immediately after dataset fetch and before polar processing

### Public Interfaces / Config
- No user-facing API change is required for the first debugging pass.
- If debugging confirms start-method sensitivity, add one explicit internal config point for multiprocessing context selection rather than relying on scattered edits or commented script variants.
- Keep `fork` out of the supported-path changes for now unless forkserver is proven unfixable.

### Test Plan
- Verify artifact staging created the three copied files in `auxiliary_tools/debug/1040-py314-hang/parallel/artifacts/`.
- Run the local forkserver repro with `polar` only and confirm the new logs show:
  - effective worker count
  - multiprocessing context
  - PID-tagged entry/exit around dataset fetch, mask fetch, and name-years attr reads
- Confirm whether the last successful log line moves from the current ambiguous region to one precise function boundary.
- Compare `num_workers=1/2/4` behavior to determine whether the failure is concurrency-sensitive or intrinsic to forkserver startup.
- If a child still exits abruptly, capture the exact last PID-tagged step before `BrokenProcessPool`.

### Assumptions
- The artifacts to stage are only the open bash file and the two open log files, not `qa.py` or `dataset_xr.py`.
- The staging destination is `auxiliary_tools/debug/1040-py314-hang/parallel/artifacts/`.
- Execution is deferred until we leave Plan Mode; this plan is the exact implementation/debugging sequence to follow next.
