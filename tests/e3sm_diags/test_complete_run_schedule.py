"""Static validation for the NERSC complete-run scheduler assets."""

from __future__ import annotations

from pathlib import Path

COMPLETE_RUN_ROOT = Path(__file__).parents[1] / "complete_run"


def test_scrontab_template_has_required_cron_controller_directives():
    template = (COMPLETE_RUN_ROOT / "complete-run.scrontab.template").read_text(
        encoding="utf-8"
    )

    for directive in (
        "--account",
        "--qos=cron",
        "--constraint=cron",
        "--time",
        "--output",
        "--open-mode=append",
    ):
        assert f"#SCRON {directive}" in template
    assert "complete-run-controller.sh" in template
    assert "0 6 * * 0" in template


def test_controller_explicitly_initializes_conda_clears_slurm_and_serializes():
    controller = (COMPLETE_RUN_ROOT / "complete-run-controller.sh").read_text(
        encoding="utf-8"
    )

    assert 'source "$CONDA_BASE/etc/profile.d/conda.sh"' in controller
    assert 'export CARTOPY_DATA_DIR="${CARTOPY_DATA_DIR:-}"' in controller
    assert 'cd "$REPOSITORY"' in controller
    assert 'conda activate "$CONTROLLER_ENV_PREFIX"' in controller
    assert "unset SLURM_MEM_PER_CPU SLURM_OPEN_MODE" in controller
    assert 'for variable in "${!SLURM_@}"' in controller
    assert "flock -n" in controller
    assert "tests.complete_run.automation" in controller
    assert "tests.complete_run.report publish" in controller
    assert "date +%V" in controller
    assert "ISO_WEEK % 2 != 0" in controller
    assert 'WORKTREE_ROOT="${WORKTREE_ROOT:-$PSCRATCH' in controller
    assert 'ENVIRONMENT_ROOT="${ENVIRONMENT_ROOT:-$PSCRATCH' in controller
    assert '"failure_count"] > 0' in controller


def test_makefile_exposes_safe_scrontab_management_commands():
    makefile = (Path(__file__).parents[2] / "Makefile").read_text(encoding="utf-8")

    for target in (
        "complete-run-scron-config",
        "complete-run-ops-init",
        "complete-run-ops-env-create",
        "complete-run-ops-env-update",
        "complete-run-ops-env-show",
        "complete-run-scron-validate",
        "complete-run-scron-install",
        "complete-run-scron-show",
        "complete-run-scron-remove",
    ):
        assert f"{target}:" in makefile
    assert "CONFIRM=YES" in makefile
