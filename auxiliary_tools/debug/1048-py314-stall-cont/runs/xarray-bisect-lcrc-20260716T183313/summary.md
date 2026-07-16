# Summary

Outcome: no stall reproduced on July 16, 2026.

Runner command:

- `bash auxiliary_tools/debug/1048-py314-stall-cont/run_xarray_bisect_qa_lcrc.sh --include-conda-releases`

Environment setup commands:

- `bash auxiliary_tools/debug/1048-py314-stall-cont/create_xarray_bisect_envs.sh`
- `bash auxiliary_tools/debug/1048-py314-stall-cont/create_xarray_latest_env.sh`

Run shape:

- Chrysalis/LCRC via `srun --pty --nodes=1 --time=01:00:00 /bin/bash`
- Per-env timeout: `30m`
- `E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1`
- `E3SM_DIAGS_REPRO_RUNS=10` was exported, but the current `parallel-lcrc/min-scripts/qa.py` still executes one diagnostics run per env

Environments tested:

- `ed_1048_xr_before_018ad08b`: xarray commit `8d271fb393372bcd2ed6ab60c9f469a1625a4aed`, Python `3.14.6`, xarray `0.0.0`
- `ed_1048_xr_after_018ad08b`: xarray commit `018ad08b12e8471b8bcc0135ce59b227f50da54b`, Python `3.14.6`, xarray `0.0.0`
- `ed_1048_xr_before_0a2d81c7`: xarray commit `43edfa34300b3513659551980a30eef393925928`, Python `3.14.6`, xarray `0.0.0`
- `ed_1048_xr_after_0a2d81c7`: xarray commit `0a2d81c7a17aab867aed362b0882d34cb89e1311`, Python `3.14.6`, xarray `0.0.0`
- `ed_1048_xr_2025120`: conda `xarray=2025.12.0`, Python `3.14.6`, xarray `2025.12.0`
- `ed_1048_xr_2026010`: conda `xarray=2026.01.0`, Python `3.14.6`, xarray `2026.1.0`
- `ed_1048_xr_latest_2026070`: conda `xarray=2026.07.0`, Python `3.14.6`, xarray `2026.7.0`

Result:

- All 7 environments finished with `status=completed`
- All 7 environments finished with `exit_code=0`
- All 7 environments had `stall_signal=no`

Conclusion: this run did not reproduce the netCDF file lock stall in either the 4 git-bisect envs or the 3 conda xarray release envs.
