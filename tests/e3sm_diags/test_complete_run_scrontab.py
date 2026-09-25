"""Tests for complete-run scrontab configuration management."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

from tests.complete_run import scrontab


def _config(tmp_path: Path) -> Path:
    config_path = tmp_path / "controller.env"
    config_path.write_text(
        "\n".join(
            (
                f"REPOSITORY={tmp_path / 'repository'}",
                f"LOG_DIR={tmp_path / 'logs'}",
                f"CONDA_BASE={tmp_path / 'conda'}",
                f"CONTROLLER_ENV_PREFIX={tmp_path / 'controller-env'}",
                f"RESULTS_ROOT={tmp_path / 'results'}",
                "SLURM_ACCOUNT=e3sm",
                "SLURM_QOS=regular",
                "SCRON_CPUS=2",
                "SCRON_MEMORY_PER_CPU=2G",
                "E3SM_DIAGS_REPOSITORY_ID=R_1",
                "E3SM_DIAGS_CATEGORY_ID=C_1",
                f"E3SM_DIAGS_TOKEN_FILE={tmp_path / 'token'}",
                "",
            )
        ),
        encoding="utf-8",
    )
    return config_path


def test_create_config_copies_template_with_private_permissions(tmp_path: Path):
    config_path = scrontab.create_config(tmp_path / "nested" / "controller.env")

    assert config_path.read_text(
        encoding="utf-8"
    ) == scrontab._CONFIG_TEMPLATE.read_text(encoding="utf-8")
    content = config_path.read_text(encoding="utf-8")
    assert "Static operations deployment location" in content
    assert "Optional local overrides" in content
    assert config_path.stat().st_mode & 0o777 == 0o600
    with pytest.raises(FileExistsError):
        scrontab.create_config(config_path)


def test_initialize_operations_clones_once_and_creates_external_config(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []
    monkeypatch.setattr(
        scrontab.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )

    checkout, config = scrontab.initialize_operations(
        tmp_path / "operations", "https://example/e3sm_diags.git", "feature-branch"
    )

    assert checkout == tmp_path / "operations" / "e3sm_diags"
    assert config.exists()
    assert (
        f"CONTROLLER_ENV_PREFIX={tmp_path / 'operations' / 'controller-env'}"
        in config.read_text(encoding="utf-8")
    )
    content = config.read_text(encoding="utf-8")
    assert (
        "REPOSITORY=/global/cfs/projectdirs/e3sm/e3sm_diags/operations/e3sm_diags"
        in content
    )
    assert "LOG_DIR=/global/cfs/projectdirs/e3sm/e3sm_diags/operations/logs" in content
    assert (tmp_path / "operations" / "logs").is_dir()
    assert calls[0][0][0][0:5] == [
        "git",
        "clone",
        "--branch",
        "feature-branch",
        "--single-branch",
    ]


def test_initialize_operations_refuses_to_replace_existing_config(tmp_path: Path):
    operations_dir = tmp_path / "operations"
    checkout = operations_dir / "e3sm_diags"
    checkout.mkdir(parents=True)
    (checkout / ".git").mkdir()
    (operations_dir / "controller.env").write_text("existing\n", encoding="utf-8")

    with pytest.raises(FileExistsError, match="controller configuration"):
        scrontab.initialize_operations(
            operations_dir, "https://example/repo.git", "main"
        )


def test_validate_config_renders_all_scheduler_placeholders(tmp_path: Path):
    config_path = _config(tmp_path)

    rendered = scrontab.validate_config(config_path)

    assert "{{" not in rendered
    assert "#SCRON --account=e3sm" in rendered
    assert "#SCRON --cpus-per-task=2" in rendered
    assert "#SCRON --mem-per-cpu=2G" in rendered
    assert str(config_path.resolve()) in rendered
    assert rendered.count("#SCRON --account=e3sm") == 4
    assert "complete-run-reporter-%j.out" in rendered
    assert "0 13 * * 0" in rendered
    assert "0 14 * * 0" in rendered
    assert "0 16 * * 1" in rendered
    assert "0 17 * * 1" in rendered


def test_validate_config_rejects_unresolved_or_relative_values(tmp_path: Path):
    config_path = _config(tmp_path)
    config_path.write_text(
        config_path.read_text(encoding="utf-8").replace(
            f"REPOSITORY={tmp_path / 'repository'}", "REPOSITORY=relative"
        ),
        encoding="utf-8",
    )

    with pytest.raises(ValueError, match="absolute path: REPOSITORY"):
        scrontab.validate_config(config_path)


def test_install_scrontab_submits_the_validated_rendering(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        if args[0] == ["scrontab", "-l"]:
            return subprocess.CompletedProcess(args[0], 0, "0 1 * * * other-job\n", "")
        return subprocess.CompletedProcess(args[0], 0, "", "")

    monkeypatch.setattr(scrontab.subprocess, "run", run)

    scrontab.install_scrontab(_config(tmp_path))

    assert calls[0][0] == (["scrontab", "-l"],)
    assert calls[1][0] == (["scrontab"],)
    assert "0 1 * * * other-job" in calls[1][1]["input"]
    assert scrontab._MANAGED_BEGIN in calls[1][1]["input"]


def test_install_replaces_only_existing_managed_block(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []
    existing = (
        "0 1 * * * other-job\n\n"
        f"{scrontab._MANAGED_BEGIN}\nold schedule\n{scrontab._MANAGED_END}\n"
    )

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        if args[0] == ["scrontab", "-l"]:
            return subprocess.CompletedProcess(args[0], 0, existing, "")
        return subprocess.CompletedProcess(args[0], 0, "", "")

    monkeypatch.setattr(scrontab.subprocess, "run", run)
    scrontab.install_scrontab(_config(tmp_path))

    installed = calls[1][1]["input"]
    assert "other-job" in installed
    assert "old schedule" not in installed
    assert installed.count(scrontab._MANAGED_BEGIN) == 1


def test_remove_scrontab_retains_unmanaged_entries(monkeypatch: pytest.MonkeyPatch):
    calls = []
    existing = (
        "0 1 * * * other-job\n\n"
        f"{scrontab._MANAGED_BEGIN}\nmanaged\n{scrontab._MANAGED_END}\n"
    )

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(args[0], 0, existing, "")

    monkeypatch.setattr(scrontab.subprocess, "run", run)
    scrontab.remove_scrontab()

    assert calls[1][0] == (["scrontab"],)
    assert calls[1][1]["input"] == "0 1 * * * other-job\n"


def test_remove_cli_invokes_managed_scrontab_removal(monkeypatch: pytest.MonkeyPatch):
    removed = []
    monkeypatch.setattr(scrontab, "remove_scrontab", lambda: removed.append(True))

    assert scrontab.main(["remove"]) == 0
    assert removed == [True]


def test_read_scrontab_accepts_common_absent_table_message(
    monkeypatch: pytest.MonkeyPatch,
):
    monkeypatch.setattr(
        scrontab.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 1, "", "no crontab for user"
        ),
    )

    assert scrontab._read_scrontab() == ""


def test_managed_scrontab_rejects_malformed_markers():
    with pytest.raises(ValueError, match="malformed"):
        scrontab._replace_managed_block(scrontab._MANAGED_BEGIN, "new")


def test_create_controller_environment_uses_configured_prefix(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []
    config_path = _config(tmp_path)

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        return subprocess.CompletedProcess(args[0], 0, "")

    monkeypatch.setattr(
        scrontab.subprocess,
        "run",
        run,
    )

    scrontab.create_controller_environment(config_path)

    assert calls[0][0][0][1:4] == ["env", "create", "--prefix"]
    assert calls[1][0][0][1:6] == [
        "run",
        "--prefix",
        str(tmp_path / "controller-env"),
        "python",
        "-m",
    ]


def test_update_controller_environment_exports_before_updating(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    calls = []
    config_path = _config(tmp_path)
    (tmp_path / "controller-env").mkdir()

    def run(*args, **kwargs):
        calls.append((args, kwargs))
        stdout = "name: controller\n" if kwargs.get("capture_output") else ""
        return subprocess.CompletedProcess(args[0], 0, stdout)

    monkeypatch.setattr(scrontab.subprocess, "run", run)

    export_path = scrontab.update_controller_environment(config_path, confirmed=True)

    assert export_path.read_text(encoding="utf-8") == "name: controller\n"
    assert calls[0][0][0][1:3] == ["env", "export"]
    assert calls[1][0][0][1:4] == ["env", "update", "--prune"]
    assert calls[-2][0][0][-2:] == ["tests.complete_run.automation", "--help"]
    assert calls[-1][0][0][-2:] == ["tests.complete_run.reporter", "--help"]


def test_update_controller_environment_requires_confirmation(tmp_path: Path):
    with pytest.raises(ValueError, match="without confirmation"):
        scrontab.update_controller_environment(_config(tmp_path), confirmed=False)


def test_update_controller_environment_refuses_an_active_controller(tmp_path: Path):
    config_path = _config(tmp_path)
    (tmp_path / "controller-env").mkdir()

    with scrontab._controller_lock(tmp_path / "results"):
        with pytest.raises(RuntimeError, match="controller or reporter is active"):
            scrontab.update_controller_environment(config_path, confirmed=True)


def test_show_controller_environment_returns_prefix_and_python_version(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    monkeypatch.setattr(
        scrontab.subprocess,
        "run",
        lambda *args, **kwargs: subprocess.CompletedProcess(
            args[0], 0, "Python 3.14.0\n"
        ),
    )

    metadata = scrontab.show_controller_environment(_config(tmp_path))

    assert metadata == {
        "prefix": str(tmp_path / "controller-env"),
        "specification": str(tmp_path / "repository" / "conda-env" / "ci.yml"),
        "python_version": "Python 3.14.0",
    }
