# July 17, 2026: Python 3.14 and Xarray Stalling Experiment

## Overview

The tests compare Python 3.14.3 and 3.14.6 across selected Xarray releases and
bisected commits to determine whether the stall is associated with a specific:

1. Python 3.14 version
2. Xarray release
3. Xarray commit
4. Another dependency (e.g., NumPy, Pandas, Dask, Matplotlib, Cartopy, xESMF, ESMPy, xgcm)

> [!NOTE]
> The NetCDF3 climatology `lock=False` workaround is disabled with
> `export E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1` so it does not mask the
> stall. See the related [commit](https://github.com/E3SM-Project/e3sm_diags/commit/5c29460ca4df2409c81bdfdc95abaf4c8d3ce75b).

## Key Findings

> [!IMPORTANT]
> A Python [change](https://docs.python.org/3/whatsnew/changelog.html) between
> Python 3.14.3 and 3.14.6 may have fixed the stall; the cause is unknown.

- **Python 3.14.3 (Min Test Run):** ❌ Stalled with Xarray >=2026.1.0
- **Python 3.14.6 (Min Test Run):** ✅ No stalls across all Xarray versions and bisected commits
- **Python 3.14.3 (Full Run):** ❌ Stalled with Xarray==2026.07.0
- **Python 3.14.6 (Full Run):** ✅ No stall with Xarray==2026.07.0
- **Test environment dependencies mostly align -- differences are limited to Python and Xarray (as expected), and a few unrelated dependencies. Rules out other dependencies as the cause of the stall.**

| Python | Xarray | Outcome |
| --- | --- | --- |
| 3.14.3 | 2025.12.0 | ✅ **Passed:** 3/3 iterations completed |
| 3.14.3 | 2026.1.0 and 2026.7.0 | ❌ **Stalled:** force-terminated after 60 minutes with exit code 143 |
| 3.14.6 | Four bisected commits and all three releases | ✅ **Passed:** 3/3 iterations completed in all seven environments |

### Recommended Action

> [!TIP]
> Keep the `lock=False` workaround so all Python 3.14 versions can run with
> Xarray >=2026.1.0 without stalling using NetCDF3 files. For read-only NetCDF3 files,
> lock=False is generally safe because they do not use HDF5,
> so xarray’s HDF5-related serialization is unnecessary and can cause avoidable contention or deadlocks.

## Test Environments Comparison

### Core Dependency Summary

**Summary:**

- **9 of 11 core dependencies match** across all environments; **Python and Xarray differ as intended**.
- `filelock`, Arrow/PyArrow, `tomlkit`, and `tzdata` differ by Python cohort, not by passing or stalling outcome.

Below are the environments used for the test cases.

| Dependency | Python 3.14.3 environments | Python 3.14.6 environments |
| --- | --- | --- |
| Python | 3.14.3 | 3.14.6 |
| Xarray releases | 2025.12.0, 2026.1.0, 2026.7.0 | 2025.12.0, 2026.1.0, 2026.7.0 |
| Xarray bisect builds | Not tested | Four Git commits installed from PyPI metadata as version 0.0.0 |
| Uxarray | 2026.06.0 | 2026.06.0 |
| NumPy | 2.4.6 | 2.4.6 |
| Pandas | 3.0.3 | 3.0.3 |
| Dask | 2026.7.1 | 2026.7.1 |
| Matplotlib | 3.11.0 | 3.11.0 |
| Cartopy | 0.24.0 | 0.24.0 |
| xESMF | 0.9.2 | 0.9.2 |
| ESMPy | 8.9.1 | 8.9.1 |
| xgcm | 0.10.0 | 0.10.0 |

### Dependency Artifacts

| CSV | Contents |
| --- | --- |
| [Core environment dependencies](core_environment_dependencies.csv) | Versions, resolved package names, and channels for 11 core packages across all 10 environments |
| [Environment dependency differences](environment_dependency_differences.csv) | The 13 packages whose version or channel differs between environments |

---

## 1. Minimum Example Test Cases

Each test runs the minimum [qa.py](auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/min-scripts/qa.py) reproduction three times per environment, with a 60-minute timeout for each run.

### Case 1: Python 3.14.3 Across Xarray Releases

**Xarray versions:** 2025.12.0, 2026.01.0, and 2026.07.0.

**Goal:** Reproduce the stalling observed with Python 3.14.3 and Xarray >=2026.01.0 (see [stalling issue context in PR #1048](https://github.com/E3SM-Project/e3sm_diags/pull/1048#issuecomment-4173567554))

**Provenance:**

| Step | Command or path |
| --- | --- |
| Create environments | `bash auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/1_create_xarray_release_envs_py3143.sh` |
| Run test matrix | `bash auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/1_run_xarray_release_qa_lcrc_py3143.sh --repro-runs 3 --timeout 60m` |
| Results | `auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/runs/xarray-release-py3143-lcrc-20260717T130132` |

**Results:**

- ✅ Xarray **2025.12.0 completed all three iterations**
- ❌ Xarray **2026.1.0 and 2026.7.0 stalled** and required force termination.
  - The stall was intermittent in onset, occurring in iteration 3 for 2026.1.0 and iteration 1 for 2026.7.0.

| Environment | Xarray source | Xarray version | Python | Result |
| --- | --- | --- | --- | --- |
| `ed_1048_xr_2025120_py3143` | Conda release | 2025.12.0 | 3.14.3 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-release-py3143-lcrc-20260717T130132/logs/ed_1048_xr_2025120_py3143.run.log) |
| `ed_1048_xr_2026010_py3143` | Conda release | 2026.1.0 | 3.14.3 | [❌ Stalled in iteration 3/3; force-terminated after 60 minutes, exit code 143](runs/xarray-release-py3143-lcrc-20260717T130132/logs/ed_1048_xr_2026010_py3143.run.log) |
| `ed_1048_xr_latest_2026070_py3143` | Conda release | 2026.7.0 | 3.14.3 | [❌ Stalled in iteration 1/3; force-terminated after 60 minutes, exit code 143](runs/xarray-release-py3143-lcrc-20260717T130132/logs/ed_1048_xr_latest_2026070_py3143.run.log) |

### Case 2: Python 3.14.6 Across Xarray Releases and Bisected Commits

**Goal:** Determine whether the stalling is fixed with a newer Python version and/or whether the stalling is associated with a specific Xarray release or commit.

**Xarray versions:** 2025.12.0, 2026.01.0, and 2026.07.0, plus four bisected commits from the Xarray 2026.01.0 development history.

**Provenance:**

| Step | Command or path |
| --- | --- |
| Create environments | `bash auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/2_create_xarray_bisect_envs_py3146.sh` |
| Run test matrix | `bash auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/2_run_xarray_bisect_qa_lcrc_py3146.sh --include-conda-releases --repro-runs 3 --timeout 60m` |
| Results | `auxiliary_tools/debug/1048-py314-stall-cont/test-scripts/runs/xarray-bisect-lcrc-20260717T132834` |

**Results:**

- ✅ All seven environments **completed** all three iterations with exit code 0 and no stall.
- Neither bisected Xarray change introduced a failure boundary under Python 3.14.6, and all three release versions passed.

| Environment | Xarray source | Xarray commit/version | Python | Result |
| --- | --- | --- | --- | --- |
| `ed_1048_xr_before_018ad08b` | Git commit | [`8d271fb393372bcd2ed6ab60c9f469a1625a4aed`](https://github.com/pydata/xarray/commit/8d271fb393372bcd2ed6ab60c9f469a1625a4aed) | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_before_018ad08b.run.log) |
| `ed_1048_xr_after_018ad08b` | Git commit | [`018ad08b12e8471b8bcc0135ce59b227f50da54b`](https://github.com/pydata/xarray/commit/018ad08b12e8471b8bcc0135ce59b227f50da54b) | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_after_018ad08b.run.log) |
| `ed_1048_xr_before_0a2d81c7` | Git commit | [`43edfa34300b3513659551980a30eef393925928`](https://github.com/pydata/xarray/commit/43edfa34300b3513659551980a30eef393925928) | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_before_0a2d81c7.run.log) |
| `ed_1048_xr_after_0a2d81c7` | Git commit | [`0a2d81c7a17aab867aed362b0882d34cb89e1311`](https://github.com/pydata/xarray/commit/0a2d81c7a17aab867aed362b0882d34cb89e1311) | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_after_0a2d81c7.run.log) |
| `ed_1048_xr_2025120` | Conda release | 2025.12.0 | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_2025120.run.log) |
| `ed_1048_xr_2026010` | Conda release | 2026.1.0 | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_2026010.run.log) |
| `ed_1048_xr_latest_2026070` | Conda release | 2026.7.0 | 3.14.6 | [✅ Completed 3/3 iterations, exit code 0, no stall](runs/xarray-bisect-lcrc-20260717T132834/logs/ed_1048_xr_latest_2026070.run.log) |

## 2. Full Run Tests Cases

These full diagnostic runs compare Python versions while holding Xarray at 2026.07.0.

### Case 1: Python 3.14.3 With Xarray 2026.07.0

**Goal:** Check if stalling occurs still with the older Python (3.14.3) and latest version of Xarray (2026.07.0)

**Provenance:**

| Item | Value |
| --- | --- |
| Command | `sbatch auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/bash-scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.py314_ed_1048_xr_latest_2026070_py3143.bash` |
| Job ID | `1253670` |
| Environment file | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_latest_2026070_py3143/model_vs_obs_1985-2014_units/prov/environment.yml` |
| Log file | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_latest_2026070_py3143/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |

**Results:**

| Status | Start | Last progress | Evidence |
| --- | --- | --- | --- |
| ⚠️ Suspected stall; job still running | 2026-07-17 15:18:12 CDT | 2026-07-17 16:20:50 CDT | At 16:54:12, the log had not changed for 33 minutes. |

The last log entry completed a GPCP 2D diagnostic. The job will be terminated if no further progress occurs.

### Case 2: Python 3.14.6 With Xarray 2026.07.0

**Goal:** Check if stalling occurs still with the latest version Python (3.14.6) and/or Xarray (2026.07.0).

**Provenance:**

| Item | Value |
| --- | --- |
| Command | `sbatch auxiliary_tools/debug/1048-py314-stall-cont/parallel-lcrc/bash-scripts/e3sm_diags_atm_monthly_180x360_aave_model_vs_obs_1985-2014.py314_ed_1048_xr_latest_2026070_py3146.bash` |
| Job ID | `1253671` |
| Environment file | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_latest_2026070_py3146/model_vs_obs_1985-2014_units/prov/environment.yml` |
| Log file | `/lcrc/group/e3sm/public_html/ac.tvo/ed_1048_xr_latest_2026070_py3146/model_vs_obs_1985-2014_units/prov/e3sm_diags_run.log` |

**Results:**

| Status | Start | Diagnostics completed | Duration |
| --- | --- | --- | --- |
| ✅ Diagnostics completed; viewer generated | 2026-07-17 15:18:12 CDT | 2026-07-17 16:01:58 CDT | About 43 minutes |

The diagnostic step completed with exit code 0. The Slurm job then exited with
code 11 because the script's post-run `rsync` used `${results_dir}` instead of
the generated `${results_dir}_units` directory. This copy failure is unrelated
to the diagnostic result.
