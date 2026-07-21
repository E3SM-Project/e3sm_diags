# July 21, 2026: Xarray Lock Stall Across Python Versions

## Overview

These tests hold Xarray at 2026.7.0 while varying Python across 3.13.12,
3.13.13, 3.13.14, and 3.14.0 through 3.14.6. They compare the E3SM Diags
NetCDF3 climatology `lock=False` workaround in two modes:

- Disabled with `E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=1`.
- Enabled with `E3SM_DIAGS_DISABLE_CLIMO_LOCK_WORKAROUND=0`.

Both modes were tested with a five-iteration minimum QA reproduction and the
full ATM monthly 180x360 model-versus-observations workflow.

## Key Findings

> [!IMPORTANT]
> With the workaround disabled, Python 3.14.1 through 3.14.4 stall in both the
> minimum and full tests. Python 3.13.12 through 3.13.14, Python 3.14.0, and
> Python 3.14.5 through 3.14.6 complete. With the workaround enabled, every
> tested environment completes both test types.

| Python | Xarray | Minimum, workaround disabled | Full, workaround disabled | Minimum, workaround enabled | Full, workaround enabled |
| --- | --- | --- | --- | --- | --- |
| 3.13.12 | 2026.7.0 | ✅ 5/5 | ✅ 00:45:18 | ✅ 5/5 | ✅ 00:47:15 |
| 3.13.13 | 2026.7.0 | ✅ 5/5 | ✅ 00:46:08 | ✅ 5/5 | ✅ 00:47:34 |
| 3.13.14 | 2026.7.0 | ✅ 5/5 | ✅ 00:45:46 | ✅ 5/5 | ✅ 00:47:11 |
| 3.14.0 | 2026.7.0 | ✅ 5/5 | ✅ 00:45:56 | ✅ 5/5 | ✅ 00:47:34 |
| 3.14.1 | 2026.7.0 | ❌ Stalled in iteration 1 | ❌ 04:00:00 timeout | ✅ 5/5 | ✅ 00:47:42 |
| 3.14.2 | 2026.7.0 | ❌ Stalled in iteration 3 | ❌ 04:00:00 timeout | ✅ 5/5 | ✅ 00:47:26 |
| 3.14.3 | 2026.7.0 | ❌ Stalled in iteration 2 | ❌ 04:00:00 timeout | ✅ 5/5 | ✅ 00:46:59 |
| 3.14.4 | 2026.7.0 | ❌ Stalled in iteration 2 | ❌ 04:00:00 timeout | ✅ 5/5 | ✅ 00:46:52 |
| 3.14.5 | 2026.7.0 | ✅ 5/5 | ✅ 00:45:56 | ✅ 5/5 | ✅ 00:47:25 |
| 3.14.6 | 2026.7.0 | ✅ 5/5 | ✅ 00:46:14 | ✅ 5/5 | ✅ 00:47:09 |

The disabled-mode failures define two Python boundaries:

1. The stall appears after Python 3.14.0 and is absent again by Python 3.14.5.
2. The enabled-mode results demonstrate that the climatology lock workaround avoids the failure
across the complete tested range, including all four affected Python releases.

## Test Environments

| Environment | Python | Xarray |
| --- | --- | --- |
| `ed_1048_xr_2026070_py31312` | 3.13.12 | 2026.7.0 |
| `ed_1048_xr_2026070_py31313` | 3.13.13 | 2026.7.0 |
| `ed_1048_xr_2026070_py31314` | 3.13.14 | 2026.7.0 |
| `ed_1048_xr_2026070_py3140` | 3.14.0 | 2026.7.0 |
| `ed_1048_xr_2026070_py3141` | 3.14.1 | 2026.7.0 |
| `ed_1048_xr_2026070_py3142` | 3.14.2 | 2026.7.0 |
| `ed_1048_xr_2026070_py3143` | 3.14.3 | 2026.7.0 |
| `ed_1048_xr_2026070_py3144` | 3.14.4 | 2026.7.0 |
| `ed_1048_xr_2026070_py3145` | 3.14.5 | 2026.7.0 |
| `ed_1048_xr_2026070_py3146` | 3.14.6 | 2026.7.0 |

## Provenance

| Test | Workaround | Results |
| --- | --- | --- |
| Five-iteration minimum QA | Disabled | [Summary](runs/python-release-xr2026070-lcrc-20260721T121830/summary.tsv) |
| Five-iteration minimum QA | Enabled | [Summary](runs/python-release-xr2026070-lcrc-20260721T123034/summary.tsv) |
| Full diagnostics | Disabled | [Summary](full-runs/python-release-xr2026070-lcrc-20260721T122041/summary.tsv) |
| Full diagnostics | Enabled | [Summary](full-runs/python-release-xr2026070-lcrc-20260721T151208/summary.tsv) |

The environments were created with
`1_create_python_release_envs_xr2026070.sh`. The minimum matrices were run with
`2_run_python_release_qa_lcrc_xr2026070.sh`, using five reproductions per
environment. The full matrices were run with
`3_run_python_release_full_lcrc_xr2026070.sh`.

## 1. Minimum Example Test Cases

Each environment runs the minimum
[qa.py](../../parallel-lcrc/min-scripts/qa.py) reproduction five times. The
jobs execute concurrently in separate LCRC Slurm allocations.

### Workaround Disabled

| Python | Job ID | Result | Slurm duration |
| --- | --- | --- | --- |
| 3.13.12 | `1255741` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py31312.run.log) | 00:44:28 |
| 3.13.13 | `1255742` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py31313.run.log) | 00:44:31 |
| 3.13.14 | `1255743` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py31314.run.log) | 00:45:05 |
| 3.14.0 | `1255744` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3140.run.log) | 00:44:40 |
| 3.14.1 | `1255745` | [❌ Stalled in iteration 1](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3141.run.log) | 01:10:01 |
| 3.14.2 | `1255746` | [❌ Completed 2/5; stalled in iteration 3](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3142.run.log) | 01:10:01 |
| 3.14.3 | `1255747` | [❌ Completed 1/5; stalled in iteration 2](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3143.run.log) | 01:10:01 |
| 3.14.4 | `1255748` | [❌ Completed 1/5; stalled in iteration 2](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3144.run.log) | 01:10:01 |
| 3.14.5 | `1255749` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3145.run.log) | 00:44:09 |
| 3.14.6 | `1255750` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T121830/logs/ed_1048_xr_2026070_py3146.run.log) | 00:44:51 |

### Workaround Enabled

All ten jobs completed all five iterations with exit code 0.

| Python | Job ID | Result | Slurm duration |
| --- | --- | --- | --- |
| 3.13.12 | `1255774` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py31312.run.log) | 00:43:19 |
| 3.13.13 | `1255776` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py31313.run.log) | 00:43:28 |
| 3.13.14 | `1255775` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py31314.run.log) | 00:43:53 |
| 3.14.0 | `1255778` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3140.run.log) | 00:43:37 |
| 3.14.1 | `1255777` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3141.run.log) | 00:42:46 |
| 3.14.2 | `1255779` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3142.run.log) | 00:43:38 |
| 3.14.3 | `1255780` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3143.run.log) | 00:43:45 |
| 3.14.4 | `1255781` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3144.run.log) | 00:43:40 |
| 3.14.5 | `1255782` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3145.run.log) | 00:44:27 |
| 3.14.6 | `1255783` | [✅ Completed 5/5](runs/python-release-xr2026070-lcrc-20260721T123034/logs/ed_1048_xr_2026070_py3146.run.log) | 00:43:49 |

The disabled-mode stall is nondeterministic in when it appears: the affected
jobs reached between zero and two completed iterations before timing out. The
workaround-enabled jobs show no corresponding failures over 50 total
iterations.

## 2. Full Run Test Cases

The full matrix runs the ATM monthly 180x360 model-versus-observations
workflow once in each environment. Each Slurm job has a four-hour walltime.

### Workaround Disabled

| Python | Job ID | Diagnostic result | Total run time |
| --- | --- | --- | --- |
| 3.13.12 | `1255762` | ✅ Viewer generated | 00:45:18 |
| 3.13.13 | `1255763` | ✅ Viewer generated | 00:46:08 |
| 3.13.14 | `1255764` | ✅ Viewer generated | 00:45:46 |
| 3.14.0 | `1255765` | ✅ Viewer generated | 00:45:56 |
| 3.14.1 | `1255766` | ❌ Stalled; Slurm timeout | 04:00:00 |
| 3.14.2 | `1255767` | ❌ Stalled; Slurm timeout | 04:00:00 |
| 3.14.3 | `1255768` | ❌ Stalled; Slurm timeout | 04:00:00 |
| 3.14.4 | `1255769` | ❌ Stalled; Slurm timeout | 04:00:00 |
| 3.14.5 | `1255770` | ✅ Viewer generated | 00:45:56 |
| 3.14.6 | `1255771` | ✅ Viewer generated | 00:46:14 |

### Workaround Enabled

| Python | Job ID | Diagnostic result | Total run time |
| --- | --- | --- | --- |
| 3.13.12 | `1255845` | ✅ Viewer generated | 00:47:15 |
| 3.13.13 | `1255846` | ✅ Viewer generated | 00:47:34 |
| 3.13.14 | `1255847` | ✅ Viewer generated | 00:47:11 |
| 3.14.0 | `1255848` | ✅ Viewer generated | 00:47:34 |
| 3.14.1 | `1255849` | ✅ Viewer generated | 00:47:42 |
| 3.14.2 | `1255850` | ✅ Viewer generated | 00:47:26 |
| 3.14.3 | `1255851` | ✅ Viewer generated | 00:46:59 |
| 3.14.4 | `1255852` | ✅ Viewer generated | 00:46:52 |
| 3.14.5 | `1255853` | ✅ Viewer generated | 00:47:25 |
| 3.14.6 | `1255854` | ✅ Viewer generated | 00:47:09 |

The six successful disabled-mode runs completed in 45-47 minutes. The four
affected releases made partial progress but produced no viewer and remained
stuck until Slurm enforced the four-hour limit. In contrast, all ten
workaround-enabled runs generated viewers in 47-48 minutes.

The completed diagnostic jobs subsequently exited with code 11 because the
inherited batch script attempts to copy `${results_dir}` instead of the
generated `${results_dir}_units` directory. This post-run `rsync` failure is
unrelated to the diagnostic outcome; viewer generation is used to classify
these runs as completed.

## 3. Python Release Interpretation

The expanded matrix narrows the affected Python interval to 3.14.1 through
3.14.4. Python 3.14.0 is an important non-stalling control, while Python
3.14.5 is the first later patch release that does not reproduce the problem.

Python 3.14.5 reverted the incremental garbage collector used in early Python
3.14 releases to the Python 3.13 generational collector. The
[Python 3.14.5 release page](https://www.python.org/downloads/release/python-3145/),
[What's New GC section](https://docs.python.org/3/whatsnew/3.14.html#garbage-collection),
and [discussion thread](https://discuss.python.org/t/reverting-the-incremental-gc-in-python-3-14-and-3-15/107014)
make that rollback a plausible explanation for why the problem disappears in
3.14.5. However, the successful Python 3.14.0 runs show that the incremental
collector alone does not explain why the stall first appears in Python 3.14.1.
The exact Python change that introduced the affected behavior remains
unidentified.

Relevant references:

- [gh-142516](https://github.com/python/cpython/issues/142516): GC issue with
  linked changes restoring the generational collector in Python 3.14.
- [gh-148144](https://github.com/python/cpython/issues/148144): incremental-GC
  frame-initialization fix.
- [Python 3.14.5rc1 changelog](https://docs.python.org/3/whatsnew/changelog.html#python-3-14-5-release-candidate-1):
  GC-related changes included before the final release.

## Conclusion

- Without the climatology lock workaround, the stall reproduces on Python
  3.14.1 through 3.14.4 and does not reproduce on Python 3.13.12 through
  3.13.14, Python 3.14.0, or Python 3.14.5 through 3.14.6.
- The same affected range appears in the five-iteration minimum matrix and the
  full diagnostics matrix.
- Enabling the workaround eliminates the observed stalls in all 20 enabled-mode
  jobs, covering 50 minimum-test iterations and 10 full diagnostic runs.
- The Python 3.14.5 GC rollback remains a plausible explanation for the upper
  boundary, but it does not by itself explain the transition between Python
  3.14.0 and 3.14.1.
