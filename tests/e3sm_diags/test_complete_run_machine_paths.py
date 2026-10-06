"""Tests for machine-selected complete-run deployment defaults."""

from __future__ import annotations

from configparser import Error as ConfigParserError
from pathlib import Path
from types import SimpleNamespace

import pytest

from tests.complete_run import machine_paths, ops, scrontab


@pytest.mark.parametrize("machine", ["pm-cpu", "pm-gpu", "unknown"])
def test_machine_mapping(machine: str, monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(
        machine_paths,
        "MachineInfo",
        lambda **kwargs: SimpleNamespace(machine=machine),
    )
    expected = machine_paths.MACHINE_PATHS.get(machine)
    assert machine_paths.detect_machine_paths() == expected
    assert (
        machine_paths.default_results_root() == machine_paths.NERSC_PATHS.results_root
    )


@pytest.mark.parametrize(
    "error", [ValueError, RuntimeError, OSError, ImportError, ConfigParserError]
)
def test_detection_failure_and_explicit_config(
    error: type[Exception], tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    def fail(**kwargs: object) -> None:
        raise error("detection unavailable")

    monkeypatch.setattr(machine_paths, "MachineInfo", fail)
    assert machine_paths.detect_machine_paths() is None
    config = tmp_path / "explicit.env"
    config.touch()
    assert ops.resolve_config(config) == config


def test_mache_unavailable(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(machine_paths, "MachineInfo", None)
    assert machine_paths.detect_machine_paths() is None


def test_machine_config_discovery_is_read_only_and_retains_precedence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkout = tmp_path / "dev" / "e3sm_diags"
    checkout.mkdir(parents=True)
    root = tmp_path / "operations"
    root.mkdir()
    config = root / "controller.env"
    config.touch()
    defaults = machine_paths.MachinePaths(root, tmp_path / "results")
    monkeypatch.setattr(ops, "_CHECKOUT", checkout)
    monkeypatch.setattr(ops, "detect_machine_paths", lambda: defaults)
    monkeypatch.delenv("E3SM_DIAGS_OPS_CONFIG", raising=False)
    before = sorted(tmp_path.rglob("*"))
    assert ops.resolve_config() == config
    assert sorted(tmp_path.rglob("*")) == before
    local = checkout.parent / "controller.env"
    local.touch()
    assert ops.resolve_config() == local
    monkeypatch.setenv("E3SM_DIAGS_OPS_CONFIG", str(config))
    assert ops.resolve_config() == config
    assert ops.resolve_config(local) == local
    with pytest.raises(FileNotFoundError, match="missing.env"):
        ops.resolve_config(tmp_path / "missing.env")
    monkeypatch.setenv("E3SM_DIAGS_OPS_CONFIG", str(tmp_path / "missing.env"))
    with pytest.raises(FileNotFoundError, match="missing.env"):
        ops.resolve_config()


@pytest.mark.parametrize("mapped", [True, False])
def test_missing_config_guidance(
    mapped: bool, tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    defaults = machine_paths.MachinePaths(tmp_path / "ops", tmp_path / "results")
    monkeypatch.setattr(ops, "_CHECKOUT", tmp_path / "dev" / "checkout")
    monkeypatch.setattr(
        ops, "detect_machine_paths", lambda: defaults if mapped else None
    )
    monkeypatch.delenv("E3SM_DIAGS_OPS_CONFIG", raising=False)
    selected = defaults.operations_dir if mapped else tmp_path / "dev"
    with pytest.raises(FileNotFoundError, match=str(selected / "controller.env")):
        ops.resolve_config()
    assert list(tmp_path.iterdir()) == []


def test_init_defaults_and_override(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    defaults = machine_paths.MachinePaths(tmp_path / "ops", tmp_path / "results")
    monkeypatch.setattr(ops, "detect_machine_paths", lambda: defaults)
    monkeypatch.delenv("OPS_INPUT_OPERATIONS_DIR", raising=False)
    calls: list[Path] = []

    def initialize(root: Path, url: str, branch: str) -> tuple[Path, Path]:
        calls.append(root)
        return root / "e3sm_diags", root / "controller.env"

    monkeypatch.setattr(scrontab, "initialize_operations", initialize)
    assert ops.main(["init"]) == 0
    override = tmp_path / "custom"
    assert ops.main(["init", "--operations-dir", str(override)]) == 0
    assert calls == [defaults.operations_dir, override]
    monkeypatch.setattr(ops, "detect_machine_paths", lambda: None)
    with pytest.raises(SystemExit):
        ops.main(["init"])
    assert ops.main(["init", "--operations-dir", str(override)]) == 0
    assert list(tmp_path.iterdir()) == []


def test_selected_results_root_in_new_configuration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    defaults = machine_paths.MachinePaths(tmp_path / "ops", tmp_path / "results space")
    monkeypatch.setattr(machine_paths, "detect_machine_paths", lambda: defaults)
    assert machine_paths.default_results_root() == defaults.results_root
    config = scrontab.create_config(tmp_path / "controller.env")
    assert scrontab._read_config(config)["RESULTS_ROOT"] == str(defaults.results_root)
    assert config.stat().st_mode & 0o777 == 0o600
    config.write_text("RESULTS_ROOT=/custom/results\n", encoding="utf-8")
    with pytest.raises(FileExistsError):
        scrontab.create_config(config)
    assert scrontab._read_config(config)["RESULTS_ROOT"] == "/custom/results"
