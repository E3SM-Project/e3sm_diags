# Plan: Automate NERSC complete-run validation (#1084)

## Scope and invariants

Automate a periodic **environment-regression** complete run from an exact
`main` revision: create a fresh CI environment, run on Perlmutter, compare with
`latest-main`, publish a deterministic report, and create one SimBoard
Discussion per comparison. The automation must retain candidate artifacts for
human review, including when diagnostics, comparison, or Slurm fail.

It must not promote a baseline, use `--allow-non-main`, loosen tolerances, or
classify a failed comparison as acceptable. Those restrictions in
`AGENTS.md`, `.claude/skills/complete-run/SKILL.md`, and the current Layer 4
documentation remain in force.

The existing configured result root is
`/global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test`; the issue's
`/global/cfs/projectdirs/...` path must not be adopted without confirmation.

## Phase 0: Confirm deployment contracts

Resolve these operational requirements before adding a scheduler or credential
integration:

1. Use NERSC's supported Slurm crontab, `scrontab`, as the scheduler. It runs
   scheduled controller jobs on NERSC's login-node pool; traditional
   `crontab` is intentionally unavailable. GitHub-hosted Actions runners are
   not needed and cannot invoke Perlmutter `sbatch` directly.
2. Confirm the canonical CFS root and its NERSC Portal URL mapping.
3. Confirm the SHA-qualified environment-name convention (for example,
   `e3sm_diags_ci_<short-main-sha>`), cleanup/retention policy, Slurm account,
   QoS, walltime, and notification owner.
4. Provision a non-personal token that can create discussions in
   `E3SM-Project/simboard`; record its storage location outside the repository.
   Obtain the GraphQL node ID for the `Complete Test Run Reports` category.
5. Decide whether all terminal outcomes are posted. The recommended policy is
   to post passed, failed, and incomplete runs so scheduler failures are
   visible; none of them trigger a promotion.

## Phase 1: Repeatable NERSC execution and local report generation

1. Add `tests/complete_run/automation.py` as an explicitly invoked,
   NERSC-login-node CLI. It will:
   - fetch and resolve `origin/main` to an immutable SHA and use a detached
     worktree, preventing a dirty checkout or a moving branch from changing a
     run;
   - create the SHA-qualified environment from that revision's
     `conda-env/ci.yml`, install that revision, and submit a CPU `sbatch` job;
   - run `tests.complete_run.run` followed by the existing comparison command
     with JSON, PNG, and HTML artifacts enabled; and
   - monitor the submitted job and preserve a machine-readable terminal status
     for submission errors, nonzero diagnostics exits, comparison failures,
     cancellations, and timeouts.

2. Add `tests/complete_run/report.py` to render deterministic JSON and Markdown
   from the run manifest, `comparison-report.json`, and the orchestration/Slurm
   status. It will include the exact result and comparison paths, their public
   URLs, Git SHA, selected sets, environment provenance, comparison status, and
   ordered failure-category counts. Missing artifacts must be reported as an
   incomplete run rather than causing the reporting step to crash.

3. Update `tests/complete_run/run.py` to export the active Conda environment to
   `prov/environment.yml` before publishing the immutable run manifest. This
   satisfies the existing documentation and lets the current comparison code
   report provenance differences reliably.

4. Update `docs/source/dev_guide/testing.rst` with the non-scheduled CLI,
   configuration inputs, artifact/report locations, and the explicit rule that
   failures are reported for human judgment only.

### Phase 1 tests and validation

- Add `tests/e3sm_diags/test_complete_run_automation.py` with mocked Git,
  Conda, Slurm, and clock interactions for command construction, immutable SHA
  selection, environment naming, submission parsing, and each terminal failure
  path.
- Add `tests/e3sm_diags/test_complete_run_report.py` for passed, comparison
  failed, diagnostics failed, timed-out, and missing-artifact report inputs;
  assert stable Markdown/JSON content and CFS-to-Portal URL conversion.
- Extend the existing complete-run run tests for environment-export success and
  failure behavior.
- Run the focused unit tests and pre-commit. After NERSC access is confirmed,
  exercise the CLI against an existing harmless result or a controlled job;
  do not run a baseline promotion.

## Phase 2: SimBoard publishing with retry safety

1. Extend `tests/complete_run/report.py` (or add a narrowly scoped publisher
   module if the renderer remains pure) with an explicit `publish` command.
   Use GitHub GraphQL's discussion-creation mutation with runtime-provided
   repository/category IDs and a token read from the approved secret location.
   Do not write tokens to CFS artifacts, logs, generated scripts, or command
   lines.

2. Persist an atomic publication receipt beside the comparison report, keyed by
   the comparison directory. A rerun must reuse the existing Discussion URL
   rather than opening duplicate threads. If posting fails, retain the rendered
   Markdown and report a publish-specific nonzero status so a later retry can
   proceed without re-running diagnostics.

3. Include the Discussion URL and publication outcome in the deterministic
   report without changing the underlying comparison verdict.

### Phase 2 tests and validation

- Mock GraphQL success, API failure, malformed response, and an existing
  receipt; verify request payloads, redaction, and idempotent retry behavior.
- Verify that a failed publish does not delete or rewrite complete-run and
  comparison artifacts.
- Use a non-production/test discussion category only if maintainers provide
  one; otherwise validate with mocked API responses.

## Phase 3: Biweekly deployment and operations

1. Add a versioned `scrontab` template and controller wrapper. Give the
   controller job the required `#SCRON` account, cron QoS, cron constraint,
   walltime, output, and append-mode settings. The wrapper must use absolute
   paths, initialize Conda explicitly, and clear inherited `SLURM_*`
   variables (including `SLURM_MEM_PER_CPU` and `SLURM_OPEN_MODE`) before it
   submits the Perlmutter diagnostics job, as required by NERSC's `scrontab`
   guidance.

2. Make the scheduler invoke the Phase 1 orchestration command and Phase 2
   publisher only after report generation. Serialize runs so overlapping
   `scrontab` invocations cannot share an environment, result directory, or
   discussion receipt. Monitor scheduled controller jobs with
   `squeue --me -q cron -O JobID,EligibleTime`.

3. Document operational ownership: manual dispatch/run command, logs, result
   retention, environment cleanup, recovery after scheduler or publication
   failure, and how maintainers review differences. State explicitly that
   baseline promotion remains a separate, explicitly confirmed manual action.

### Phase 3 validation

- Validate the cron/workflow syntax and shell wrapper with dry-run or mocked
  commands.
- Perform one maintainer-approved end-to-end scheduled execution on NERSC,
  confirm its public report and exactly one Discussion thread, then inspect the
  artifacts. Stop and report if comparison fails; do not promote.

## Affected files

- New: `tests/complete_run/automation.py`
- New: `tests/complete_run/report.py`
- Modified: `tests/complete_run/run.py`
- New: `tests/e3sm_diags/test_complete_run_automation.py`
- New: `tests/e3sm_diags/test_complete_run_report.py`
- Modified: relevant existing `tests/e3sm_diags/test_complete_run_*.py`
- Modified: `docs/source/dev_guide/testing.rst`
- Phase 3 only: NERSC `scrontab` template and controller wrapper

## Phase 4: Persistent operations controller environment

### Goal

Move the login-node-only controller environment out of home storage and into
the persistent operations directory. Keep it deliberately maintained rather
than updating it from `scrontab`. The fresh diagnostic environment remains a
separate timestamped, SHA-qualified `$PSCRATCH` environment solved from the
exact `origin/main` revision's `conda-env/ci.yml` on every run.

### Configuration and controller changes

1. Replace the name-only `CONTROLLER_ENV` configuration contract with an
   absolute `CONTROLLER_ENV_PREFIX`, normally
   `/global/cfs/projectdirs/e3sm/e3sm_diags/operations/controller-env`.
   Validate that the prefix is absolute and use `conda activate <prefix>` in
   the controller wrapper.
2. Update `complete-run-operations-init` so it writes that prefix into the
   external mode-0600 controller configuration. It must continue to refuse to
   replace an existing checkout or configuration file.
3. Document the separation of concerns:
   - the persistent CFS controller prefix runs only login-node orchestration;
   - every diagnostics job gets a newly solved `$PSCRATCH` environment from the
     exact main revision's `ci.yml`;
   - the CFS run manifest and `prov/environment.yml`, not the PSCRATCH prefix,
     are durable diagnostic-environment provenance.

### Explicit controller-environment lifecycle commands

1. Add `make complete-run-controller-env-create CONFIG=<controller.env>`.
   It reads the external configuration, refuses an existing prefix, creates it
   from the operations checkout's `conda-env/ci.yml`, and installs that checkout.
2. Add
   `make complete-run-controller-env-update CONFIG=<controller.env> CONFIRM=YES`.
   This command is manual and must never be invoked by `scrontab`. It must:
   - require `CONFIRM=YES`;
   - acquire the same controller lock used by scheduled runs, refusing an
     update while a controller is active;
   - export the existing controller prefix to a timestamped CFS provenance file
     before mutation;
   - run `conda env update --prune` using the operations checkout's `ci.yml`;
   - reinstall the checkout; and
   - verify `python -m tests.complete_run.automation --help` through the
     updated prefix.
3. Add `make complete-run-controller-env-show CONFIG=<controller.env>` to
   display the configured prefix and its Conda metadata without modification.
4. Keep removal out of automatic operations. If a removal command is added,
   require `CONFIRM=YES` and preserve/export provenance first.

### Tests and validation

- Add unit tests for absolute-prefix validation, bootstrap rendering, create
  refusal on an existing prefix, update command construction, lock contention,
  environment export ordering, and post-update verification failure.
- Add Makefile static/dry-run tests for create, show, and confirmation-gated
  update targets.
- Update scheduler/controller tests to assert prefix activation rather than a
  named environment.
- Run focused complete-run scheduler tests, `bash -n` on the controller,
  Make dry-runs, documentation build, and `pre-commit run --all-files`.
- After review, have an operations owner create the controller prefix, record
  its export, and perform one approved controller-environment update while no
  `scrontab` controller is active. Do not run a baseline promotion.
