---
name: complete-run
description: Run and interpret E3SM Diags Layer 4 complete-run validation on NERSC Perlmutter. Use when validating a high-risk change or a release against the accepted complete-run baseline, when checking for dependency regressions with a freshly solved conda environment, or when re-comparing or promoting complete-run results.
---

# Layer 4 Complete-Run Validation

Runs a large cross-section of diagnostics against HPC-hosted data and compares
the resulting netCDF files and PNG plots with the accepted `latest-main`
baseline.

**The procedure lives in `docs/source/dev_guide/testing.rst`**, sections "Layer
4: Complete-Run Validation" and "Complete-Run Validation": environment
construction, the sbatch submission, the Makefile targets, and the tolerances.
Follow it there rather than reproducing it here. This file covers the judgment
the procedure does not: which mode to run, what not to do, and how to read the
result.

This is a manual, human-in-the-loop workflow. It is expensive and it writes to
shared space on CFS.

## Guardrails

1. **Never run `make promote-complete` or `baseline promote` without explicit
   confirmation from the user in the current conversation.** Promotion
   overwrites `latest-main`, the baseline every other developer compares
   against. A prior "yes" for a different run does not carry over.
2. **Never use `--allow-non-main`.** It exists for a maintainer overriding the
   branch check by hand. Surface it, do not reach for it.
3. **Stop and report after a failed comparison.** Show the user the summary and
   the diff artifacts. Do not decide whether differences are acceptable, do not
   loosen tolerances to make a comparison pass, and do not promote.
4. **Prefer re-comparing over re-running.** The diagnostics take most of an
   hour; the comparison takes minutes and can be repeated against the same
   results directory as often as needed.

## Pick a mode, and say which one you picked

The comparison cannot tell you whether a numerical difference came from the
code or from the dependencies. Decide which one you are holding fixed.

**Mode A — code validation.** Validating a branch. Rebuild the baseline's
environment from its `prov/environment.yml` so any difference is attributable
to code. The default for PR validation.

**Mode B — environment regression.** Checking whether current conda-forge
breaks the plots. Clean `main` checkout, fresh solve from `conda-env/dev.yml`.
`dev.yml` carries floating constraints that only a fresh solve exercises. The
pre-release gate, not a per-PR step.

Mode B assumes the code matches the baseline's. `main` moves, so confirm that
before trusting the attribution:

```bash
python -m tests.complete_run.baseline show --channel main   # note the manifest sha
git diff --stat <baseline-sha>..HEAD -- e3sm_diags/
```

An empty diff means the shipped package is identical and differences really are
dependencies. A non-empty one means you are running Mode A whether you meant to
or not; say so.

## Preflight

```bash
git status --porcelain -uno          # expect clean
git rev-parse --abbrev-ref HEAD      # Mode B requires main
python -m tests.complete_run.baseline show --channel main
```

`baseline show` prints the resolved baseline directory and its manifest. The
manifest's `environment` block records the python version, platform, and
curated package versions the baseline was produced with; that is what your run
will be diffed against. If it reports a missing or broken link, nothing has
been promoted yet — stop and tell the user the comparison has nothing to run
against.

`run.py` validates input paths before starting, so a bad path fails in seconds
rather than an hour into the run.

## Read the report

At `<results-parent>/comparison/<dev>-vs-<baseline>/comparison-report.json`, in
this order:

1. **`environment`** — first, before anything else. In Mode B this section *is*
   the result. In Mode A it should be empty; if it is not, the environment did
   not reproduce the baseline's and numerical differences are not attributable
   to code. `conda_environment` always differs when the run used a fresh
   environment, and that difference alone is not a finding.
2. **`status`** and **`summary.failure_count`** — the verdict.
3. **The per-category lists** — `missing_dev_files`, `missing_baseline_files`,
   `missing_variables`, `nan_location_mismatches`, `shape_mismatches`,
   `tolerance_failures`, `missing_dev_images`, `missing_baseline_images`,
   `image_mismatches`.
4. **`paths.diff_artifact_dir`** — the rendered PNG diffs.

**Then stop.** Summarize the failure categories, counts, and the most
significant diff artifacts, and let the user judge. Do not characterize
differences as expected or acceptable on your own.

## Handling a failure

A failed comparison leaves the candidate results and artifacts in place. Re-run
the comparison, not the diagnostics — see `testing.rst` for the targets.

Varying tolerances is a review aid, not a way to pass. Do it only at the user's
request, and when reporting a comparison that passed at a loosened tolerance,
say that is what happened.

**If the job hit the walltime**, the results directory has no
`baseline_manifest.json` — it is written after the diagnostics complete, so its
absence means the run was cut short. Such a directory can still be compared,
but its environment provenance degrades to "unavailable" and it can never be
promoted. Report it and offer to remove it rather than comparing partial
output; a rerun goes to a new timestamped directory regardless.

## Promotion — only on explicit confirmation

Preconditions: the changes are reviewed, approved, and merged to `main`; the
run came from a `main` checkout; and **the user has just told you, in this
conversation, to promote.**

If the run is not on `main` the promotion is refused. That refusal is correct.
Report it and stop.

## Judgment calls the CLI will not make for you

- **A `--set` subset is not a validation.** Comparing a subset run against a
  full baseline reports every unrun set as missing files. Useful for iterating
  on one diagnostic, useless as a gate.
- **`--results-dir` is used verbatim**, skipping the timestamp-branch-commit
  suffix that the default path in `params.py` gets. A second run into the same
  path fails on the immutable-manifest guard. Pass a path not used before.
- **Do not edit `params.py` to change a run.** Those defaults are the shared
  contract that makes results comparable. Use the `run.py` CLI flags;
  `python -m tests.complete_run.run --help` lists them.
