# July 20, 2026: Xarray Lock Stall Across Python Versions

## Overview

These tests hold Xarray at 2026.7.0 while varying Python across 3.13.14,
3.14.3, 3.14.4, 3.14.5, and 3.14.6. The goal is to identify which Python
patch releases reproduce the E3SM Diags stall first observed with Python
3.14.3 and Xarray >=2026.1.0.

Both the minimum QA matrix and the full diagnostic runs use the same five
Conda environments.

> [!NOTE]
> The NetCDF3 climatology `lock=False` workaround is disabled with
> `export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1` so it does not mask the
> stall.

## Key Findings

> [!IMPORTANT]
> The minimum test identifies a boundary between Python 3.14.4 and 3.14.5.
> Python 3.14.3 and 3.14.4 stalled, while Python 3.13.14, 3.14.5, and 3.14.6
> completed all three iterations. The most plausible Python-side change is the
> Python 3.14.5 garbage-collector reversion: the incremental GC shipped in
> Python 3.14.0-3.14.4 was reverted to the Python 3.13 generational GC. Refer to [Python 3.14.5 Changelog Review](#3-python-3145-changelog-review) for details.

| Test | Python | Xarray | Outcome |
| --- | --- | --- | --- |
| Minimum | 3.13.14 | 2026.7.0 | ✅ Passed all 3 iterations |
| Minimum | 3.14.3 | 2026.7.0 | ❌ Stalled in iteration 1; timed out after 60 minutes |
| Minimum | 3.14.4 | 2026.7.0 | ❌ Passed iteration 1, stalled in iteration 2; timed out after 60 minutes |
| Minimum | 3.14.5 | 2026.7.0 | ✅ Passed all 3 iterations |
| Minimum | 3.14.6 | 2026.7.0 | ✅ Passed all 3 iterations |
| Full | 3.13.14 | 2026.7.0 | ✅ Diagnostics completed; viewer generated in 44:54 |
| Full | 3.14.3 | 2026.7.0 | ❌ Confirmed stalled; last progress is an unmatched `open_mfdataset start` |
| Full | 3.14.4 | 2026.7.0 | ❌ Confirmed stalled; last progress is an unmatched `open_mfdataset start` |
| Full | 3.14.5 | 2026.7.0 | ✅ Diagnostics completed; viewer generated in 44:01 |
| Full | 3.14.6 | 2026.7.0 | ✅ Diagnostics completed; viewer generated in 44:58 |

The full-run results are consistent with the minimum test. Python 3.13.14,
3.14.5, and 3.14.6 completed in approximately 44-45 minutes. Python 3.14.3 and
3.14.4 remained active in Slurm at **2026-07-20 15:38 CDT**, but their
provenance logs had not advanced since **14:52:58 CDT** and **14:47:21 CDT**,
respectively. These jobs were manually terminated with `scancel 1255367 1255368`.

## Test Environments

| Environment | Python | Xarray |
| --- | --- | --- |
| `ed_1048_xr_2026070_py31314` | 3.13.14 | 2026.7.0 |
| `ed_1048_xr_2026070_py3143` | 3.14.3 | 2026.7.0 |
| `ed_1048_xr_2026070_py3144` | 3.14.4 | 2026.7.0 |
| `ed_1048_xr_2026070_py3145` | 3.14.5 | 2026.7.0 |
| `ed_1048_xr_2026070_py3146` | 3.14.6 | 2026.7.0 |

## Provenance

| Step | Command or path |
| --- | --- |
| Create environments | `bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/1_create_python_release_envs_xr2026070.sh` |
| Run minimum matrix | `bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/2_run_python_release_qa_lcrc_xr2026070.sh --repro-runs 3 --timeout 60m` |
| Run full matrix | `bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/3_run_python_release_full_lcrc_xr2026070.sh` |
| Python 3.14.5 changelog review | https://www.python.org/downloads/release/python-3145/, https://docs.python.org/3/whatsnew/3.14.html#garbage-collection, https://docs.python.org/3/whatsnew/changelog.html#python-3-14-5-release-candidate-1, https://docs.python.org/3/whatsnew/changelog.html#python-3-14-5-final, and https://discuss.python.org/t/reverting-the-incremental-gc-in-python-3-14-and-3-15/107014 |
| GC/RSS/stall confirmation runner | `bash auxiliary_tools/debug/1048-py314-stall-cont/investigations/260720-xarray-lock-stall-across-py-versions/4_run_py3144_gc_trace_min_lcrc.sh` |
| Minimum results | [Minimum-run summary](runs/python-release-xr2026070-lcrc-20260720T133718/summary.tsv) |
| Full results | [Full-run summary](full-runs/python-release-xr2026070-lcrc-20260720T133755/summary.tsv) |

## 1. Minimum Example Test Cases

Each environment runs the minimum
[qa.py](../../parallel-lcrc/min-scripts/qa.py) reproduction three times. The
entire three-iteration command has a 60-minute timeout. Runs execute
concurrently in separate LCRC Slurm allocations.

### Results

| Environment | Python | Job ID | Result | Slurm duration |
| --- | --- | --- | --- | --- |
| `ed_1048_xr_2026070_py31314` | 3.13.14 | `1255361` | [✅ Completed 3/3 iterations, exit code 0](runs/python-release-xr2026070-lcrc-20260720T133718/logs/ed_1048_xr_2026070_py31314.run.log) | 00:26:26 |
| `ed_1048_xr_2026070_py3143` | 3.14.3 | `1255362` | [❌ Stalled in iteration 1, timeout exit code 124](runs/python-release-xr2026070-lcrc-20260720T133718/logs/ed_1048_xr_2026070_py3143.run.log) | 01:00:17 |
| `ed_1048_xr_2026070_py3144` | 3.14.4 | `1255363` | [❌ Completed iteration 1, stalled in iteration 2, timeout exit code 124](runs/python-release-xr2026070-lcrc-20260720T133718/logs/ed_1048_xr_2026070_py3144.run.log) | 01:00:17 |
| `ed_1048_xr_2026070_py3145` | 3.14.5 | `1255364` | [✅ Completed 3/3 iterations, exit code 0](runs/python-release-xr2026070-lcrc-20260720T133718/logs/ed_1048_xr_2026070_py3145.run.log) | 00:26:42 |
| `ed_1048_xr_2026070_py3146` | 3.14.6 | `1255365` | [✅ Completed 3/3 iterations, exit code 0](runs/python-release-xr2026070-lcrc-20260720T133718/logs/ed_1048_xr_2026070_py3146.run.log) | 00:26:33 |

The successful iterations each took approximately 8.5-8.8 minutes:

| Python | Iteration 1 | Iteration 2 | Iteration 3 |
| --- | --- | --- | --- |
| 3.13.14 | 530.2 s | 512.7 s | 507.3 s |
| 3.14.3 | ❌ Stalled | ❌ Not run | ❌ Not run |
| 3.14.4 | 521.8 s | ❌ Stalled | ❌ Not run |
| 3.14.5 | 522.5 s | 516.9 s | 527.4 s |
| 3.14.6 | 519.1 s | 516.9 s | 521.0 s |

## 2. Full Run Test Cases

The full matrix runs the ATM monthly 180x360 model-versus-observations
workflow once in each environment. Each Slurm job has a four-hour walltime.

> [!WARNING]
> The Python 3.14.3 and 3.14.4 jobs were manually cancelled after they were
> confirmed stalled. When checked at 2026-07-20 15:38 CDT, their provenance logs
> showed no progress for more than 45 minutes after the last unmatched
> `open_mfdataset start`.

### Results as of 2026-07-20 15:38 CDT

| Environment | Python | Job ID | Diagnostic result | Total run time | Last provenance event |
| --- | --- | --- | --- | --- | --- |
| `ed_1048_xr_2026070_py31314` | 3.13.14 | `1255366` | ✅ Viewer generated | 00:44:54 | `Viewer HTML generated` at 14:23:01 CDT |
| `ed_1048_xr_2026070_py3143` | 3.14.3 | `1255367` | ❌ Confirmed stalled; manually cancelled | 01:34:05 at 15:38 CDT before cancellation | `open_mfdataset start` for `LWCF` at 14:52:58 CDT; no matching `done` |
| `ed_1048_xr_2026070_py3144` | 3.14.4 | `1255368` | ❌ Confirmed stalled; manually cancelled | 01:33:58 at 15:38 CDT before cancellation | `open_mfdataset start` for `SST` at 14:47:21 CDT; no matching `done` |
| `ed_1048_xr_2026070_py3145` | 3.14.5 | `1255369` | ✅ Viewer generated | 00:44:01 | `Viewer HTML generated` at 14:49:11 CDT |
| `ed_1048_xr_2026070_py3146` | 3.14.6 | `1255370` | ✅ Viewer generated | 00:44:58 | `Viewer HTML generated` at 14:53:51 CDT |

The Python 3.14.3 and 3.14.4 provenance logs both end immediately after
`Climo backend open_mfdataset start`, with no matching `done` event. At
2026-07-20 15:38 CDT, `squeue` still reported both jobs as `RUNNING`, and the
log file modification times were unchanged at the final provenance timestamps.
Both jobs were then manually cancelled with `scancel 1255367 1255368`. The logs
also confirm that the lock workaround was disabled by the environment variable.
This is consistent with the failure mode seen in the minimum test.

### Full-Run Artifacts

| Python | Results directory | Provenance log |
| --- | --- | --- |
| 3.13.14 | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py31314/model_vs_obs_1985-2014_units` | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py31314/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |
| 3.14.3 | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3143/model_vs_obs_1985-2014_units` | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3143/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |
| 3.14.4 | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3144/model_vs_obs_1985-2014_units` | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3144/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |
| 3.14.5 | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3145/model_vs_obs_1985-2014_units` | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3145/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |
| 3.14.6 | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3146/model_vs_obs_1985-2014_units` | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_2026070_py3146/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |

The three completed diagnostics generated their viewers successfully. Their
Slurm jobs subsequently exited with code 11 because the inherited batch
script attempts to copy `${results_dir}` instead of the generated
`${results_dir}_units` directory. This post-run `rsync` failure is unrelated
to the diagnostic outcome.

## 3. Python 3.14.5 Changelog Review

Most likely cause: Python 3.14.5 reverted the incremental GC used in Python
3.14.0-3.14.4 back to the Python 3.13 generational GC. The
[Python 3.14.5 release page](https://www.python.org/downloads/release/python-3145/),
[What's New GC section](https://docs.python.org/3/whatsnew/3.14.html#garbage-collection),
and [discussion thread](https://discuss.python.org/t/reverting-the-incremental-gc-in-python-3-14-and-3-15/107014)
all cite production memory-pressure reports as the reason for the rollback.

Why it fits: with Xarray, inputs, and the lock workaround held constant, Python
3.14.3 and 3.14.4 stall; Python 3.13.14, 3.14.5, and 3.14.6 complete.

Relevant links:

- [gh-142516](https://github.com/python/cpython/issues/142516): GC issue with
  linked PRs to restore the generational GC in Python 3.14.
- [gh-148144](https://github.com/python/cpython/issues/148144): lower-confidence
  incremental-GC frame initialization fix.
- [Python 3.14.5rc1 changelog](https://docs.python.org/3/whatsnew/changelog.html#python-3-14-5-release-candidate-1):
  contains the GC-related changes.
- [Python 3.14.5 final changelog](https://docs.python.org/3/whatsnew/changelog.html#python-3-14-5-final):
  no similarly strong Linux default-build `open_mfdataset` match.

Best next confirmation: rerun Python 3.14.4 with worker-level stack dumps plus
RSS/GC telemetry. A lock wait should show the stuck worker's wait site; memory
pressure or GC pathology should show up in the RSS/GC samples while the final
event remains an unmatched `Climo backend open_mfdataset start`.

## Conclusion

- The July 20 tests narrow the Python transition associated with the stall to
**after Python 3.14.4 and by Python 3.14.5**.
- Python 3.13.14 does not reproduce the problem, indicating that it is not a general behavior of all supported
Python versions.
- The Python 3.14.5 changelog points to the GC rollback as the
most plausible fix.
- The full runs support the same boundary: Python 3.14.3 and
3.14.4 were confirmed stalled in `open_mfdataset` and manually cancelled, while
Python 3.13.14, 3.14.5, and 3.14.6 completed.
