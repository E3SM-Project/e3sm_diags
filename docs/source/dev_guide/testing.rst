Testing E3SM Diagnostics
========================

Testing Architecture
--------------------

E3SM Diagnostics uses four test layers across local, CI/CD, and LCRC
environments:

1. Unit tests for code correctness
2. Targeted image-regression tests for rendered output
3. Broad integration tests for diagnostic workflows
4. Complete-run validation against HPC-hosted data

.. figure:: _static/testing-architecture.svg
   :alt: Testing architecture diagram showing test layers by environment.

Recommended Contributor Workflow
--------------------------------

For most changes:

1. Run Layer 1 during local development.
2. Run Layer 2 for changes that may affect plots or rendered output.
3. Run Layer 3 when broader workflow coverage is needed.
4. Run the repository's default local checks before opening a pull request.
5. Let CI/CD enforce Layers 1 through 3 on the pull request.
6. Run Layer 4 manually for high-risk changes requiring full NERSC validation.

Local Test Layers
-----------------

Before running a test or baseline-promotion command, activate the E3SM Diags
Conda environment:

.. code-block:: bash

   conda activate <e3sm_diags_env>

Layer 1: Unit Tests
~~~~~~~~~~~~~~~~~~~

Unit tests check code correctness and API stability. Run them first during
local development:

.. code-block:: bash

   make test-unit

Layer 2: Targeted Image-Regression Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These tests detect pixel-level plot regressions using targeted synthetic cases
and committed PNG baselines. Run them after Layer 1 when a code or dependency
change may affect rendered output:

.. code-block:: bash

   make test-image-regression

The suite currently covers:

- ``lat_lon``
- ``polar``
- ``zonal_mean_2d``
- ``cosp_histogram``

Baselines and their dependency metadata are stored in
``tests/integration/baselines/``.

Investigating Failures
^^^^^^^^^^^^^^^^^^^^^^

Rerun the test with a persistent artifact directory:

.. code-block:: bash

   IMAGE_REGRESSION_ARTIFACT_DIR=tests/integration/image_check_failures \
       make test-image-regression

Inspect ``tests/integration/image_check_failures`` to determine whether the
change is expected. Each failed case includes:

- The generated image
- ``runtime_metadata.json``
- ``dependency_diff.json``, comparing the runtime environment with the
  committed ``baseline_metadata.json``

.. note::

   GitHub Actions uploads these artifacts when an image-regression test fails.
   Download them from the workflow run summary page.

Updating Baselines
^^^^^^^^^^^^^^^^^^

Update baselines only when the plot change is intentional.

The preferred method is the manual ``Update Image Baselines`` GitHub Actions
workflow. It regenerates baselines on ``main`` using the authoritative
``conda-env/ci.yml`` and Python 3.14 environment.

After the workflow completes, rebase affected branches onto ``main`` and rerun
CI.

To refresh baselines locally:

.. code-block:: bash

   conda env create -f conda-env/ci.yml
   conda activate e3sm_diags_ci
   python -m tests.integration.refresh_plot_image_baselines
   make test-image-regression

To refresh one case:

.. code-block:: bash

   python -m tests.integration.refresh_plot_image_baselines --case polar

Use the same ``conda-env/ci.yml`` and Python 3.14 environment as the main CI
visual-regression gate. Commit the updated PNGs and
``baseline_metadata.json``.

Layer 3: Broad Downloaded-Data Integration Tests
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

These tests exercise broader diagnostic workflows with downloaded data. They
verify that workflows complete and generate outputs, but do not perform
pixel-level image comparisons.

Run them when broader integration coverage is needed:

.. code-block:: bash

   make test-integration

The integration target downloads its required data automatically. It runs with
``CHECK_IMAGES=False``, making Layer 3 a workflow smoke test rather than the
visual-regression authority.

By default, ``tests.integration.download_data`` uses the local
``/e3sm_diags_downloaded_data`` directory when available. Otherwise, it uses
``crane export`` to copy that directory from the OCI image used by CI.

Use ``--source-root`` or ``--image`` for nonstandard data sources.

Layer 4: Complete-Run Validation
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Layer 4 runs a broad set of diagnostics against HPC-hosted data and compares
netCDF outputs and PNG plots with an accepted baseline. Use it for high-risk
changes, release validation, and scheduled environment regression checks.
See `Complete-Run Validation`_ for instructions.

Complete-Run Validation
-----------------------

Manual Validation
~~~~~~~~~~~~~~~~~

1. **Choose and create the environment on a NERSC login node.**

   Use scratch for temporary environments and CFS for durable results:

   .. code-block:: bash

      STAMP=$(date -u +%Y%m%d)-$(git rev-parse --short HEAD)
      ENV_PREFIX="$SCRATCH/e3sm_diags_complete_run_$STAMP"

   For **code validation**, reproduce the baseline's dependencies:

   .. code-block:: bash

      mamba env create -f <baseline-dir>/prov/environment.yml -p "$ENV_PREFIX"

   For **environment regression**, solve fresh dependencies:

   .. code-block:: bash

      mamba env create -f conda-env/dev.yml -p "$ENV_PREFIX"

   To isolate dependency effects, use the baseline's exact code revision.
   A newer ``main`` revision may introduce code changes as well.

   Activate the environment and install the package from the checkout:

   .. code-block:: bash

      conda activate "$ENV_PREFIX"
      pip install .

2. **Submit validation from the repository root.**

   Adjust the account, QoS, and walltime for your allocation:

   .. code-block:: bash

      cat > "$SCRATCH/complete_run_$STAMP.sbatch" <<EOF
      #!/bin/bash
      #SBATCH --account=e3sm
      #SBATCH --qos=regular
      #SBATCH --constraint=cpu
      #SBATCH --nodes=1
      #SBATCH --time=01:00:00
      #SBATCH --output=$SCRATCH/complete_run_$STAMP.log

      set -eo pipefail
      source "$(conda info --base)/etc/profile.d/conda.sh"
      conda activate "$ENV_PREFIX"
      cd "$(pwd)"
      make test-complete-validate
      EOF

      sbatch "$SCRATCH/complete_run_$STAMP.sbatch"

   Source ``conda.sh`` to enable activation in the batch shell. Run from the
   repository root to test the working tree. Request a full CPU node: the
   diagnostics default to 24 workers, and ``enso_diags`` can have memory spikes.

   Alternatively, request an interactive allocation, then run validation from
   the repository root:

   .. code-block:: bash

      salloc --nodes 1 --qos interactive --time 04:00:00 --constraint cpu --account=e3sm
      conda activate "$ENV_PREFIX"
      make test-complete-validate

3. **Review the results.**

   Each run creates an immutable timestamped directory with a branch and
   commit suffix under:

   .. code-block:: text

      /global/cfs/cdirs/e3sm/www/e3sm_diags/complete-run-test/

   Validation compares against ``latest-main`` and preserves candidate outputs,
   a JSON report, and PNG diffs. Each result also records its environment in
   ``prov/environment.yml`` and its manifest, so scratch environments can be
   removed after review.

   Check the report's ``environment`` section before interpreting differences.
   For code validation, unexpected package or platform differences prevent
   attributing output changes solely to code. For environment regression,
   review both dependency changes and their effects on outputs.

Review and Baseline Management
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Repeat or Customize a Comparison
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Repeat comparisons without rerunning diagnostics. Failed comparisons preserve
candidate results and review artifacts.

.. code-block:: bash

   make test-complete-compare RUN_DIR=<results-dir>

Override the default ``latest-main`` baseline when needed:

.. code-block:: bash

   make test-complete-compare \
       RUN_DIR=<results-dir> \
       BASELINE_DIR=<baseline-dir>

Reports and PNG diffs are saved under ``comparison/`` beside the result
directories. Comparison behavior:

* **netCDF:** relative tolerance ``1e-5`` and absolute tolerance ``0.0``.
  Override with ``--rtol`` and ``--atol`` on ``tests.complete_run.compare``.
  Disclose any pass that requires wider tolerances.
* **PNG:** images are classified as identical, cosmetic, or reviewable.
  Small rendering shifts may be cosmetic; layout, plotted-content, and compact
  text changes remain reviewable. Review newly cosmetic results before
  accepting them in an environment regression.
* **HTML:** ``--write-diff-html`` creates ``index.html`` with baseline,
  candidate, and diff images, ranked from ``STRUCTURAL`` to ``MINOR``.
  It includes diagnostic-set and severity filters, reports identical and
  cosmetic counts, and implies ``--write-diff-pngs``.

To compare only PNGs:

.. code-block:: bash

   python -m tests.complete_run.compare \
       --dev-dir <results-dir> \
       --mode images

The default pixel mismatch threshold is ``0.0002``, matching the targeted
image-regression suite. Change ``--image-mismatch-threshold`` only for a
reviewed environment difference.

Promote an Approved Baseline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After approval and merge, generate and review results from ``main`` before
promoting them:

.. code-block:: bash

   make test-complete
   make test-complete-compare RUN_DIR=<main-results-dir>
   make promote-complete RUN_DIR=<main-results-dir>

Automated Environment Regression
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The NERSC login-node controller validates an exact ``origin/main`` revision
using a detached worktree and a fresh timestamped, SHA-qualified environment
from that revision's ``ci.yml``. The environment is created in the CPU Slurm
allocation, not the memory-constrained cron controller. It preserves results,
Slurm status, JSON/PNG/HTML comparisons, and
``automation-report`` files.

The compute allocation must be able to access the configured Conda channels
and package index, or have the required packages available in its Conda cache.

The operations owner configures the account, QoS, walltime, CFS-to-Portal
mapping, retention, and notifications. Use ``$PSCRATCH`` for temporary
worktrees and diagnostic environments; retain results on CFS for review.

Run Once
^^^^^^^^

.. code-block:: bash

   python -m tests.complete_run.automation \
       --worktree-root "$PSCRATCH/e3sm_diags-worktrees" \
       --environment-root "$PSCRATCH/e3sm_diags-environments" \
       --account e3sm

Schedule Biweekly Runs
^^^^^^^^^^^^^^^^^^^^^^

Scheduled runs use the controller wrapper and
``tests/complete_run/complete-run.scrontab.template``.
The controller starts at 06:00 Pacific every Sunday and runs only on even ISO
weeks, giving a biweekly Sunday cadence. The time avoids the Sunday 02:00
daylight-saving transition and peak weekday use.

1. **Initialize the operations directory.**

   From an existing checkout:

   .. code-block:: bash

      OPS_DIR=/global/cfs/projectdirs/e3sm/e3sm_diags/operations
      make complete-run-ops-init OPERATIONS_DIR="$OPS_DIR"

   This creates a non-public operations layout with a clean controller
   checkout (``e3sm_diags/``), external configuration (``controller.env``), and
   logs (``logs/``). It defaults to ``main``, clones only if the checkout is
   absent, and never overwrites an existing configuration. To test unmerged
   automation changes, add ``BRANCH=devops/1084-automate-complete-test``.

   Keep candidate results, detached worktrees, and diagnostic environments
   outside the controller checkout.

2. **Create the E3SM Diags reporting token.**

   A repository administrator must enable Discussions in
   ``E3SM-Project/e3sm_diags`` and create the ``Complete Test Run Reports``
   category. Use a dedicated machine account with organization membership
   and write access to the repository.

   While signed in as that account, create a fine-grained token at
   https://github.com/settings/personal-access-tokens/new. Select resource
   owner ``E3SM-Project``, restrict repository access to ``e3sm_diags``, and grant
   **Discussions: read and write**. Set the expiration to the maximum allowed value 
   (366 days).

   Store the token under the account that runs the controller:

   .. code-block:: bash

      make complete-run-ops-token-create

   The command creates
   ``$HOME/.config/e3sm_diags/e3sm_diags-token`` with mode ``0600``. To use a
   different location, set ``TOKEN_FILE``:

   .. code-block:: bash

      make complete-run-ops-token-create TOKEN_FILE=/absolute/path/to/token

   Keep the token outside the repository; never put its value in
   ``controller.env``.

3. **Configure the controller.**

   .. code-block:: bash

      $EDITOR "$OPS_DIR/controller.env"

   Set the operational values, including:

   * ``CONDA_BASE``: the base installation of Conda, typically ``/global/homes/v/<user>/miniforge3``.
   * ``CONTROLLER_ENV_PREFIX``: the persistent login-node environment,
     typically ``$OPS_DIR/controller-env``. Each diagnostics job uses a
     separate fresh environment in ``$PSCRATCH``.
   * ``E3SM_DIAGS_TOKEN_FILE``: ``$HOME/.config/e3sm_diags/e3sm_diags-token``.
   * ``SLURM_ACCOUNT``, ``SLURM_QOS``, ``SLURM_NODES``,
     ``SLURM_WALLTIME``, and ``SLURM_CONSTRAINT``: diagnostics batch-job
     resources. The template configures ``e3sm``, ``regular``, ``1``,
     ``02:00:00``, and ``cpu``, respectively. Slurm selects the CPU partition
     for the regular QoS.
   * ``SCRON_CPUS`` and ``SCRON_MEMORY_PER_CPU``: resources for the login-node
     controller allocation. The cron partition permits two CPUs and at most
     ``2G`` per CPU, so the template defaults to a 4G allocation. These values
     are independent of diagnostics-job resources.

   Keep the configuration outside the repository with mode ``0600``.
   If an existing operations directory lacks configuration, create it first:

   .. code-block:: bash

      make complete-run-scron-config CONFIG="$OPS_DIR/controller.env"

4. **Create the controller environment, validate, and install the schedule.**

   Run each command only after the previous one succeeds:

   .. code-block:: bash

      cd "$OPS_DIR/e3sm_diags"
      make complete-run-ops-env-create CONFIG="$OPS_DIR/controller.env"
      make complete-run-scron-validate CONFIG="$OPS_DIR/controller.env"
      make complete-run-scron-install CONFIG="$OPS_DIR/controller.env"

Maintain Scheduled Runs
^^^^^^^^^^^^^^^^^^^^^^^

From the controller checkout, inspect the schedule and jobs:

.. code-block:: bash

   make complete-run-scron-show
   squeue --me -q cron -O JobID,EligibleTime

To remove the schedule:

.. code-block:: bash

   make complete-run-scron-remove CONFIRM=YES

Update the persistent controller environment manually, never from ``scrontab``:

.. code-block:: bash

   OPS_DIR=/global/cfs/projectdirs/e3sm/e3sm_diags/operations
   make complete-run-ops-env-show CONFIG="$OPS_DIR/controller.env"
   make complete-run-ops-env-update CONFIG="$OPS_DIR/controller.env" CONFIRM=YES

The update holds the controller lock, exports the current environment to
``operations/provenance/``, updates from the checkout's ``ci.yml``, reinstalls
the checkout, and verifies the controller CLI.

The operations owner reviews differences, manages result and environment
retention, and retries failed publication using preserved Markdown reports.
Temporary environments can also expire under the normal ``$PSCRATCH`` purge
policy.

.. important::

   Automation publishes only reviewable comparison failures to E3SM Diags
   Discussions. Clean comparisons produce no Discussion, even when environment
   provenance differs. It never promotes baselines, passes
   ``--allow-non-main``, or reinterprets failed comparisons. Baseline promotion
   requires separate, explicitly confirmed manual review and action.


CI/CD Workflows
---------------

Main CI/CD Workflow
~~~~~~~~~~~~~~~~~~~

The main GitHub Actions workflow runs on pull requests and ``main``. It
includes:

1. Layer 1 unit tests
2. Layer 2 targeted image-regression tests
3. Layer 3 integration smoke tests with ``CHECK_IMAGES=False``

Layer 2 is the authoritative visual-regression gate. Layer 3 provides broader
workflow coverage without image matching.

Manual Image-Baseline Refresh
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The manual ``Update Image Baselines`` workflow updates committed Layer 2
baselines directly on ``main``. It is intended for legitimate plot changes
caused by dependency updates or other approved changes, avoiding a
baseline-only pull request.

The workflow:

1. Regenerates all Layer 2 baselines.
2. Reruns the targeted image-regression suite.
3. Pushes to ``main`` only when the diff is limited to
   ``tests/integration/baselines/``.

Use this workflow only for intentional baseline updates. Normal code changes
must use the standard pull request workflow.

.. note::

   If verification fails, the workflow uploads the generated images, runtime
   metadata, and dependency differences for review.

E3SM-Unified Advisory Compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The manual ``E3SM Unified Latest Release Advisory Compatibility`` workflow
runs Layer 2 against the latest released Linux ``nompi`` ``e3sm-unified``
package from conda-forge.

This is an advisory production-compatibility check, not the authoritative
visual-regression gate. It compares plots from the released E3SM-Unified
environment with baselines generated by the main CI environment. Dependency
differences, such as different Matplotlib versions, may therefore cause image
mismatches without indicating an ``e3sm_diags`` regression.

Run this workflow when evaluating compatibility with the released
E3SM-Unified environment. Review uploaded artifacts before classifying a
failure as a code regression.

Implementation Details
^^^^^^^^^^^^^^^^^^^^^^

The workflow:

1. Starts from ``conda-env/ci.yml``.
2. Resolves the latest released ``e3sm-unified`` package metadata from
   ``conda-forge/linux-64/repodata.json.bz2``.
3. Substitutes the released package dependencies into the CI environment.
4. Caches Conda packages using the generated environment hash.
5. Runs the Layer 2 image-regression suite.

Baseline updates remain governed by the main Layer 2 environment on ``main``.
