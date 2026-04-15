## Plan: Diagnose Py3.14 + Xarray 2026.1 Stall

Most likely, you are hitting a close-lock lifecycle interaction, not just a generic open/read race.
Working hypothesis: with Xarray >=2026.1.0, close now takes a stricter netCDF4 lock path, and in your current flow the same dataset ownership can be closed more than once across nested helpers while running under forked multiprocessing. That can convert a previously harmless pattern into a stall.

This plan includes explicit upstream evidence collection from external PR pages before narrowing root-cause conclusions.

**Steps**
1. Map the full open -> load -> close lifecycle for climo and time-series paths in [e3sm_diags/driver/utils/dataset_xr.py](e3sm_diags/driver/utils/dataset_xr.py).
2. Confirm ownership boundaries and identify duplicate close opportunities between helper and caller paths in [e3sm_diags/driver/utils/dataset_xr.py](e3sm_diags/driver/utils/dataset_xr.py).
3. Execute upstream evidence fetch checklist (command-ready):
4. Primary fetch for xarray PR 10788:
	- Tool: fetch_webpage
	- Input URL: https://github.com/pydata/xarray/pull/10788
	- Query target: "netCDF4 close locking change, files touched, rationale"
	- Expected artifact: concise summary + list of touched backend files/functions.
5. Primary fetch for e3sm_diags commit b1e98504:
	- Tool: fetch_webpage
	- Input URL: https://github.com/E3SM-Project/e3sm_diags/commit/b1e98504f03955606339a6d373ad976cf034c1c9
	- Query target: "dataset close semantics introduced, functions changed"
	- Expected artifact: concise summary + changed close/load call sites.
6. Fallback A if fetch_webpage extraction fails:
	- Command: `curl -L https://github.com/pydata/xarray/pull/10788.patch | sed -n '1,260p'`
	- Command: `curl -L https://github.com/E3SM-Project/e3sm_diags/commit/b1e98504f03955606339a6d373ad976cf034c1c9.patch | sed -n '1,260p'`
	- Expected artifact: patch hunks with file/function anchors.
7. Fallback B for richer metadata if needed:
	- Command: `curl -s https://api.github.com/repos/pydata/xarray/pulls/10788/files`
	- Command: `curl -s https://api.github.com/repos/E3SM-Project/e3sm_diags/commits/b1e98504f03955606339a6d373ad976cf034c1c9`
	- Expected artifact: machine-readable file lists and patch snippets.
8. Record evidence in a short note with exact mappings:
	- For each upstream hunk, capture: file, function, lock/close behavior, and predicted impact on e3sm_diags close ownership.
9. Correlate upstream changes with process execution in [e3sm_diags/e3sm_diags_driver.py](e3sm_diags/e3sm_diags_driver.py) and parallel controls in [e3sm_diags/parameter/core_parameter.py](e3sm_diags/parameter/core_parameter.py).
10. Reproduce with a minimal matrix using [auxiliary_tools/debug/1048-py314-stall-cont/parallel/min-scripts/qa.py](auxiliary_tools/debug/1048-py314-stall-cont/parallel/min-scripts/qa.py):
11. Python 3.14 + Xarray 2025.12.0 vs 2026.1.0.
12. num_workers = 1 vs >1.
13. explicit close path enabled vs disabled.
14. RAW_NETCDF4_CLIMO_OPEN enabled vs disabled.
15. Determine if the stall point is specifically in close-finally after inner helper close, or in sibling-process close contention.
16. Select mitigation with your preferred direction: single-owner close contract first, version gating only if strictly needed.
17. Validate no leak/race regressions with existing close-behavior tests in [tests/e3sm_diags/driver/utils/test_dataset_xr.py](tests/e3sm_diags/driver/utils/test_dataset_xr.py).

**Relevant files**
- [e3sm_diags/driver/utils/dataset_xr.py](e3sm_diags/driver/utils/dataset_xr.py): lifecycle and close ownership.
- [e3sm_diags/e3sm_diags_driver.py](e3sm_diags/e3sm_diags_driver.py): Dask multiprocessing configuration.
- [e3sm_diags/parameter/core_parameter.py](e3sm_diags/parameter/core_parameter.py): worker and multiprocessing settings.
- [auxiliary_tools/debug/1040-py314-hang/parallel/notes.md](auxiliary_tools/debug/1040-py314-hang/parallel/notes.md): prior stall evidence/workarounds.
- [auxiliary_tools/debug/1048-py314-stall-cont/parallel/min-scripts/qa.py](auxiliary_tools/debug/1048-py314-stall-cont/parallel/min-scripts/qa.py): reproducer focus.
- [tests/e3sm_diags/driver/utils/test_dataset_xr.py](tests/e3sm_diags/driver/utils/test_dataset_xr.py): behavior contracts around close.

**Verification**
1. Capture exact stall location (open, load, or close) across the matrix.
2. Confirm upstream evidence alignment by mapping xarray PR/commit hunks to observed stall points.
3. Confirm whether removing duplicate close path eliminates stall in local experiment.
4. Re-run close-related unit tests to ensure no ownership/leak regressions.
5. Run one representative parallel diagnostics workflow for regression confidence.

**Decisions captured**
- Priority: root-cause confirmation first.
- Preferred mitigation: single-owner close contract.