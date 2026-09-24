"""Static validation for the NERSC complete-run scheduler assets."""

from __future__ import annotations

import os
import subprocess
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
    assert (
        "SLURM_ACCOUNT|SLURM_CONSTRAINT|SLURM_NODES|SLURM_QOS|SLURM_WALLTIME"
        in controller
    )
    assert "flock -n" in controller
    assert "tests.complete_run.automation" in controller
    assert "tests.complete_run.report publish" in controller
    assert "date +%V" in controller
    assert "ISO_WEEK % 2 != 0" in controller
    assert 'WORKTREE_ROOT="${WORKTREE_ROOT:-$PSCRATCH' in controller
    assert 'ENVIRONMENT_ROOT="${ENVIRONMENT_ROOT:-$PSCRATCH' in controller
    assert '"failure_count"] > 0' in controller
    assert 'if [[ ! -s "$COMPLETION_FILE" ]]' in controller
    assert "Automation did not write completion metadata" in controller
    assert "SLURM_PARTITION" not in controller
    assert '--nodes "${SLURM_NODES:-1}"' in controller
    assert '--constraint "${SLURM_CONSTRAINT:-cpu}"' in controller
    assert '--walltime "${SLURM_WALLTIME:-02:00:00}"' in controller


def test_controller_handles_missing_automation_completion_file(tmp_path: Path):
    conda_base = tmp_path / "conda"
    hook = conda_base / "etc" / "profile.d" / "conda.sh"
    hook.parent.mkdir(parents=True)
    hook.write_text("conda() { :; }\n", encoding="utf-8")
    executable_dir = tmp_path / "bin"
    executable_dir.mkdir()
    (executable_dir / "date").write_text(
        "#!/bin/sh\nprintf '02\\n'\n", encoding="utf-8"
    )
    (executable_dir / "python").write_text("#!/bin/sh\nexit 1\n", encoding="utf-8")
    for executable in executable_dir.iterdir():
        executable.chmod(0o755)

    config = tmp_path / "controller.env"
    config.write_text(
        "\n".join(
            (
                f"REPOSITORY={tmp_path}",
                f"CONDA_BASE={conda_base}",
                f"CONTROLLER_ENV_PREFIX={tmp_path / 'controller-env'}",
                f"RESULTS_ROOT={tmp_path / 'results'}",
                f"PSCRATCH={tmp_path / 'pscratch'}",
                "SLURM_ACCOUNT=e3sm",
                "E3SM_DIAGS_REPOSITORY_ID=repository",
                "E3SM_DIAGS_CATEGORY_ID=category",
                f"E3SM_DIAGS_TOKEN_FILE={tmp_path / 'token'}",
                "",
            )
        ),
        encoding="utf-8",
    )
    environment = {**os.environ, "PATH": f"{executable_dir}:{os.environ['PATH']}"}

    completed = subprocess.run(
        ["bash", str(COMPLETE_RUN_ROOT / "complete-run-controller.sh"), str(config)],
        capture_output=True,
        text=True,
        env=environment,
        check=False,
    )

    assert completed.returncode == 1
    assert "Automation did not write completion metadata" in completed.stderr


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
