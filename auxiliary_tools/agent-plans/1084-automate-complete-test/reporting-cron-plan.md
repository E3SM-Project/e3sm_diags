# Complete-Run Reporting Cron Plan

## Goal

Stop holding a cron/login allocation while a regular CPU diagnostics job waits
in the Slurm queue.

## Architecture

Split the current controller into two independent roles.

### Submission controller

Run on the existing biweekly schedule. It must:

1. Resolve `origin/main`.
2. Create the detached worktree.
3. Create an immutable automation run directory and initial `status.json`.
4. Submit the regular CPU diagnostics job.
5. Record the Slurm job ID and set `stage` to `submitted`.
6. Exit without polling, rendering reports, or publishing Discussions.

### Reporting controller

Run frequently, such as every 15 minutes. It must:

1. Inspect submitted and finalized run directories under
   `<RESULTS_ROOT>/automation/`.
2. Leave jobs still visible in `squeue` unchanged.
3. Read a batch-written final `status.json` after a job leaves `squeue`.
4. Use `sacct` only when a job left no batch-written status.
5. Render automation reports and conditionally publish comparison failures.
6. Exit promptly.

## Status Lifecycle

Use these stages:

```text
submission_failed
submitted
environment_failed
diagnostics_failed
comparison_failed
passed
cancelled
timed_out
slurm_failed
job_completed_without_status
```

`submitted` means Slurm accepted the job. It must not be reported as a
submission failure.

## Idempotency and Locking

- Retain the submission lock.
- Add a distinct reporter lock at
  `<RESULTS_ROOT>/automation/reporter.lock`.
- Skip runs with an existing final automation report.
- Reuse the existing publication receipt so repeated reporter runs never create
  duplicate Discussions.
- Retry a later reporting invocation when Slurm accounting data is not yet
  available.

## Scheduler Assets

- Keep the existing biweekly submission cron entry.
- Add a frequent reporter cron entry, for example:

  ```cron
  */15 * * * * <repository>/tests/complete_run/complete-run-reporter.sh <config>
  ```

- The reporter needs a short cron allocation because it does not create Conda
  environments, run diagnostics, compare results, or poll during queue waits.

## Implementation

1. Add a reporting-only Python entry point, preferably a dedicated
   `tests.complete_run.reporter` module.
2. Make `complete-run-controller.sh` submission-only.
3. Add `complete-run-reporter.sh` for reporting and publication.
4. Move current post-submission polling, report rendering, and publication from
   the submission wrapper to the reporter.
5. Update the controller environment template, scrontab template, and testing
   documentation.

## Test Coverage

Add tests that verify:

1. Submission exits immediately after successful `sbatch`.
2. The reporter skips jobs still queued or running.
3. A completed batch-written status produces a report.
4. A completed job without status is classified with `sacct`.
5. Accounting lag is retried rather than misclassified.
6. Repeated reporter invocations do not duplicate reports or Discussions.
7. Submission and reporter locks serialize their independent work.
8. The shell wrappers invoke their respective Python entry points.

## Operational Validation

1. Submit a one-off controller test.
2. Confirm the submission cron job exits after recording `submitted`.
3. Confirm the regular diagnostics job remains queued independently.
4. Run the reporter manually after diagnostics finish.
5. Confirm report generation and conditional publication.
6. Enable the recurring reporter cron schedule.
