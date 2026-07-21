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

- **Performance implications:** In comparable NetCDF3 validation runs, the
  workaround-enabled jobs averaged 1:28 (3.2%) slower. This small observed
  performance hit is outweighed by reliably completing affected runs in about
  47 minutes instead of timing out after four hours; separate cohorts mean it
  is not established workaround overhead.

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
| Enabled-versus-disabled output validation | Python 3.13.14 and 3.14.6 | [Validation report](validation/netcdf3-climo-lock-workaround-20260721/summary.md) |

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

### Performance Implications

Across the six versions that completed in both modes, enabled averaged 00:47:21
versus 00:45:53 disabled: 00:01:28 (3.2%) slower. For Python 3.14.1-3.14.4,
enabled averaged 00:47:15 instead of timing out at 04:00:00. Because the
cohorts ran separately, the 3.2% difference does not establish overhead.

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
The tested Conda interpreters are standard GIL builds
(`Py_GIL_DISABLED=0`), so changelog entries explicitly limited to
free-threaded builds do not apply.

The [Python 3.14 changelog](https://docs.python.org/3/whatsnew/changelog.html)
contains no entry that directly identifies an Xarray, NetCDF, file-lock, or
standard-build thread deadlock spanning exactly 3.14.1 through 3.14.4. The
following entries are the plausible candidates after filtering by release
boundary and runtime mechanism:

| Candidate | Release and possible role | Assessment |
| --- | --- | --- |
| [gh-139951](https://github.com/python/cpython/issues/139951) | 3.14.1; possible introduction-side change | Changed incremental-GC accounting to count actually tracked objects and avoid excessive collections on growing heaps. It is the most relevant standard-build runtime change at the 3.14.0/3.14.1 boundary because it changes collection cadence and therefore object-finalization timing. However, it was intended to fix a performance regression already present in 3.14.0, and the changelog reports no deadlock. **Plausibility: medium, not proven.** |
| [gh-142516](https://github.com/python/cpython/issues/142516) | 3.14.5rc1; possible fix | Replaced the incremental collector in the default build with the Python 3.13 generational collector. This exactly matches the upper boundary, applies to these GIL-enabled environments, and is a broad enough runtime change to remove a timing-sensitive lock/finalizer interaction. **Plausibility: high for the fix boundary.** |
| [gh-148144](https://github.com/python/cpython/issues/148144) | 3.14.4; incremental-GC defect fix | Initialized `_PyInterpreterFrame.visited` so incremental GC would not read an uninitialized byte from copied generator and frame objects. This confirms a real incremental-GC correctness defect, but 3.14.4 still stalls, so this change did not fix the observed issue. **Plausibility: low as the direct cause; useful supporting evidence.** |
| [gh-131788](https://github.com/python/cpython/issues/131788) | 3.14.1; possible introduction-side multiprocessing change | Made `multiprocessing` resource-tracker sends re-entrant safe. It fits the lower version boundary and changes synchronization behavior, but the reported problem concerns resource cleanup rather than mid-diagnostic climatology reads. **Plausibility: low.** |
| [gh-146313](https://github.com/python/cpython/issues/146313) | 3.14.5rc1; possible fix | Fixed a resource-tracker deadlock in `os.waitpid()` during interpreter shutdown when a forked child retained the tracker pipe. The observed E3SM Diags stalls occur during active diagnostics, not shutdown, and are eliminated by changing only the Xarray lock argument. **Plausibility: low.** |

Several entries can be deprioritized or used to rule out alternatives:

- [gh-142206](https://github.com/python/cpython/issues/142206) restored the
  Python 3.14.0 multiprocessing resource-tracker protocol in Python 3.14.2.
  The stall persists through 3.14.2-3.14.4, so the protocol introduced in
  3.14.1 is not sufficient to explain the issue.
- [gh-137017](https://github.com/python/cpython/issues/137017) changed
  `threading.Thread.is_alive()` in 3.14.1 to remain true until OS-thread cleanup
  completes. This could expose a pre-existing wait, but it does not explain why
  passing `lock=False` prevents the stall.
- Free-threaded-only fixes such as `gh-142048`, `gh-144513`, and `gh-148820`
  are excluded because the tested interpreters have the GIL enabled.

The leading hypothesis is therefore a timing-sensitive interaction between
the incremental GC and Xarray/backend lock or resource finalization. The
3.14.1 GC-accounting change may expose that interaction, and replacing the
collector in 3.14.5 may remove it. This remains an inference from matching
release boundaries, not a demonstrated CPython root cause.

The most discriminating follow-up tests are:

1. Run Python 3.14.4 with the workaround disabled but `gc.disable()` active.
   Completion would implicate cyclic-GC activity, though not a specific patch.
2. Build Python 3.14.4 with the generational-GC change from `gh-142516`
   backported and rerun the minimum test. This isolates the strongest proposed
   3.14.5 fix from all other 3.14.5 changes.
3. Bisect CPython between 3.14.0 and 3.14.1, starting with the `gh-139951`
   backport commits. This tests the only plausible introduction-side GC change
   identified in the changelog.
4. Retain worker stack dumps and GC callbacks during these runs to distinguish
   a backend lock wait from a long GC pause or finalizer-induced deadlock.

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
