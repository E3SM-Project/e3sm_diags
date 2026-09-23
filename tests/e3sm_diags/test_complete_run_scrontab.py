"""Tests for complete-run scrontab configuration management."""

from __future__ import annotations

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
                "CONTROLLER_ENV=ed_dev_1084",
                f"RESULTS_ROOT={tmp_path / 'results'}",
                "SLURM_ACCOUNT=e3sm",
                "SIMBOARD_REPOSITORY_ID=R_1",
                "SIMBOARD_CATEGORY_ID=C_1",
                f"SIMBOARD_TOKEN_FILE={tmp_path / 'token'}",
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
    assert config_path.stat().st_mode & 0o777 == 0o600
    with pytest.raises(FileExistsError):
        scrontab.create_config(config_path)


def test_validate_config_renders_all_scheduler_placeholders(tmp_path: Path):
    config_path = _config(tmp_path)

    rendered = scrontab.validate_config(config_path)

    assert "{{" not in rendered
    assert "#SCRON --account=e3sm" in rendered
    assert str(config_path.resolve()) in rendered


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
    monkeypatch.setattr(
        scrontab.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )

    scrontab.install_scrontab(_config(tmp_path))

    assert calls[0][0] == (["scrontab"],)
    assert "#SCRON --account=e3sm" in calls[0][1]["input"]
