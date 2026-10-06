"""Tests for the consolidated complete-run operator interface."""

from __future__ import annotations

import json
import os
import shlex
import subprocess
from pathlib import Path

import pytest

from tests.complete_run import ops, scrontab
from tests.e3sm_diags.test_complete_run_scrontab import _config


def test_configuration_precedence(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    checkout = tmp_path / "checkout"
    checkout.mkdir()
    monkeypatch.setattr(ops, "_CHECKOUT", checkout)
    default = _config(tmp_path)
    other = tmp_path / "other.env"
    other.write_text("", encoding="utf-8")
    monkeypatch.delenv("E3SM_DIAGS_OPS_CONFIG", raising=False)
    assert ops.resolve_config() == default
    monkeypatch.setenv("E3SM_DIAGS_OPS_CONFIG", str(other))
    assert ops.resolve_config() == other
    assert ops.resolve_config(default) == default
    with pytest.raises(FileNotFoundError, match="ops-init"):
        ops.resolve_config(tmp_path / "missing.env")


def test_help_needs_no_configuration(capsys: pytest.CaptureFixture[str]) -> None:
    assert ops.main(["help"]) == 0
    output = capsys.readouterr().out
    assert output.index("Everyday") < output.index("Administration")
    assert "ops-env ACTION=create" in output
    for alias in (
        "ops-env-create",
        "ops-env-update",
        "ops-env-show",
        "ops-schedule-show",
    ):
        assert alias not in output


@pytest.mark.parametrize(
    "command,extra",
    [
        ("run", []),
        ("report", []),
        ("enable", []),
        ("disable", []),
        ("env", ["--action", "update"]),
    ],
)
def test_confirmation_precedes_discovery(
    command: str,
    extra: list[str],
    capsys: pytest.CaptureFixture[str],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.delenv("CONFIRM", raising=False)
    with pytest.raises(SystemExit) as error:
        ops.main([command, *extra])
    assert error.value.code == 2
    assert "CONFIRM=YES" in capsys.readouterr().err


def test_environment_action_is_explicit(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    with pytest.raises(SystemExit):
        ops.main(["env", "--config", str(_config(tmp_path))])
    assert "ACTION=create or ACTION=update" in capsys.readouterr().err


def test_logs_filter_and_bound(
    tmp_path: Path, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)
    root = tmp_path / "logs"
    root.mkdir()
    (root / "complete-run-controller-123.out").write_text(
        "old\nsecond\nlast\n", encoding="utf-8"
    )
    (root / "complete-run-controller-456.out").write_text(
        "wrong job\n", encoding="utf-8"
    )
    ops.logs(config, job="123", lines=2)
    output = capsys.readouterr().out
    assert "second\nlast\n" in output
    assert "old\n" not in output
    assert "wrong job" not in output
    assert "No reporter logs found for JOB=123" in output
    with pytest.raises(ValueError):
        ops.logs(config, job="../*", lines=2)
    with pytest.raises(ValueError):
        ops.logs(config, job=None, lines=0)


def test_dashboard_is_read_only_with_missing_tools_and_malformed_metadata(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)
    run = tmp_path / "results" / "automation" / "abc123-20261001-010000"
    run.mkdir(parents=True)
    (run / "status.json").write_text("{oops", encoding="utf-8")
    (run / "automation-report.json").write_text("[]", encoding="utf-8")
    before = sorted(str(path) for path in tmp_path.rglob("*"))

    def missing(*args: object, **kwargs: object) -> None:
        raise FileNotFoundError("tool not installed")

    monkeypatch.setattr(ops.subprocess, "run", missing)
    ops.dashboard(config)
    output = capsys.readouterr().out
    assert "squeue is not installed" in output
    assert "Malformed or unreadable metadata" in output
    assert "Future recurring occurrences, not previous outcomes" in output
    assert str(config) in output
    assert sorted(str(path) for path in tmp_path.rglob("*")) == before


def test_latest_run_uses_timestamp_not_sha(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    root = tmp_path / "automation"
    for name in ("zzz-20260901-010000", "aaa-20261001-010000"):
        run = root / name
        run.mkdir(parents=True)
        (run / "status.json").write_text(
            json.dumps({"job_id": "123", "stage": "passed"}), encoding="utf-8"
        )
    monkeypatch.setattr(ops, "_inspect_command", lambda arguments: "123 COMPLETED 0:0")
    ops._latest_run({"RESULTS_ROOT": str(tmp_path)})
    output = capsys.readouterr().out
    assert "aaa-20261001" in output
    assert "zzz-20260901" not in output
    assert "123 COMPLETED" in output


def test_shortcut_quotes_paths_and_forwards_arguments(tmp_path: Path) -> None:
    root = tmp_path / "with spaces and 'quote"
    root.mkdir()
    config = _config(root)
    text = ops.shortcut(config)
    # Mock make with a shell function's external-command equivalent, not the real target.
    executable = root / "bin"
    executable.mkdir()
    fake_make = executable / "make"
    fake_make.write_text("#!/bin/bash\nprintf '<%s>\\n' \"$@\"\n", encoding="utf-8")
    fake_make.chmod(0o755)
    result = subprocess.run(
        ["bash", "-c", text + '\ne3sm-ops run CONFIRM=YES "LINES=two words"'],
        capture_output=True,
        text=True,
        check=True,
        cwd=tmp_path,
        env={**os.environ, "PATH": f"{executable}:{os.environ['PATH']}"},
    )
    assert f"<{root / 'repository'}>" in result.stdout
    assert "<ops-run>" in result.stdout
    assert "<CONFIRM=YES>" in result.stdout
    assert "<LINES=two words>" in result.stdout
    assert f"<CONFIG={config}>" in result.stdout


@pytest.mark.parametrize(
    "command,component", [("run", "controller"), ("report", "reporter")]
)
def test_manual_commands_invoke_separate_wrappers(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, command: str, component: str
) -> None:
    calls = []
    monkeypatch.setattr(
        ops.subprocess, "run", lambda *args, **kwargs: calls.append((args, kwargs))
    )
    config = _config(tmp_path)
    ops.main([command, "--config", str(config), "--confirm", "YES"])
    assert calls[0][0][0] == [
        str(
            tmp_path
            / "repository"
            / "tests"
            / "complete_run"
            / f"complete-run-{component}.sh"
        ),
        str(config),
        "--manual",
    ]


def test_enable_validates_before_install_and_disable_only_removes_managed_schedule(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []
    config = _config(tmp_path)
    monkeypatch.setattr(
        ops, "validate_deployment", lambda path: calls.append("validate")
    )
    monkeypatch.setattr(
        scrontab, "install_scrontab", lambda path: calls.append("install")
    )
    monkeypatch.setattr(scrontab, "remove_scrontab", lambda: calls.append("remove"))
    ops.main(["enable", "--config", str(config), "--confirm", "YES"])
    ops.main(["disable", "--config", str(config), "--confirm", "YES"])
    assert calls == ["validate", "install", "remove"]


def test_enable_rejects_nonexecutable_scripts(tmp_path: Path) -> None:
    config = _config(tmp_path)
    with pytest.raises(FileNotFoundError, match="executable"):
        ops.validate_deployment(config)


def test_update_uses_shared_lock(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path)
    calls = []
    monkeypatch.setattr(ops, "_fast_forward", lambda path: calls.append("update"))
    monkeypatch.setattr(
        ops, "validate_deployment", lambda path: calls.append("validate")
    )
    with scrontab._controller_lock(tmp_path / "results"):
        with pytest.raises(RuntimeError, match="active"):
            ops.update_checkout(config)
    assert not calls
    ops.update_checkout(config)
    assert calls == ["update", "validate"]


@pytest.mark.parametrize(
    "dirty,branch,ahead,expected",
    [
        (" M file", "main", "0", "dirty"),
        ("", "", "0", "detached"),
        ("", "main", "1", "ahead"),
    ],
)
def test_update_rejects_unsafe_checkout(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    dirty: str,
    branch: str,
    ahead: str,
    expected: str,
) -> None:
    def git(arguments: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        command = arguments[3:]
        if command[0] == "status":
            output = dirty
        elif command[0] == "branch":
            output = branch
        elif command[0] == "rev-list":
            output = f"{ahead} 1"
        elif "--git-path" in command:
            output = str(tmp_path / "absent")
        else:
            output = "origin/main"
        return subprocess.CompletedProcess(arguments, 0, output, "")

    monkeypatch.setattr(ops.subprocess, "run", git)
    with pytest.raises(ValueError, match=expected):
        ops._fast_forward(tmp_path)


def test_update_is_fast_forward_only_and_does_not_touch_conda(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    calls = []

    def git(arguments: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        command = arguments[3:]
        calls.append(command)
        if command[0] == "status":
            output = ""
        elif command[0] == "rev-list":
            output = "0 1"
        elif "--git-path" in command:
            output = str(tmp_path / "absent")
        elif command[0] == "config":
            output = "deployment-remote"
        else:
            output = "origin/main"
        return subprocess.CompletedProcess(arguments, 0, output, "")

    monkeypatch.setattr(ops.subprocess, "run", git)
    ops._fast_forward(tmp_path)
    assert ["fetch", "deployment-remote"] in calls
    assert ["merge", "--ff-only", "origin/main"] in calls


@pytest.mark.parametrize("component", ["controller", "reporter"])
def test_wrappers_manual_bypasses_only_schedule_guards(
    tmp_path: Path, component: str
) -> None:
    executable = tmp_path / "bin"
    executable.mkdir()
    marker = tmp_path / "invocation"
    for name, body in {
        "date": 'case "$1" in +%V) printf "01\\n" ;; *) printf "00\\n" ;; esac',
        "python": f'[[ "$ACTIVATED" == yes ]] || exit 1; printf "%s\\n" "$@" > {shlex.quote(str(marker))}',
    }.items():
        path = executable / name
        path.write_text("#!/bin/bash\n" + body + "\n", encoding="utf-8")
        path.chmod(0o755)
    conda = tmp_path / "conda" / "etc" / "profile.d"
    conda.mkdir(parents=True)
    (conda / "conda.sh").write_text(
        "conda() { export ACTIVATED=yes; }\n", encoding="utf-8"
    )
    config = _config(tmp_path)
    (tmp_path / "repository").mkdir()
    script = ops._CHECKOUT / "tests" / "complete_run" / f"complete-run-{component}.sh"
    environment = {
        **os.environ,
        "PATH": f"{executable}:{os.environ['PATH']}",
        "PSCRATCH": str(tmp_path),
    }
    subprocess.run(
        [str(script), str(config)], env=environment, check=True, capture_output=True
    )
    assert not marker.exists()
    subprocess.run(
        [str(script), str(config), "--manual"],
        env=environment,
        check=True,
        capture_output=True,
    )
    assert (
        "tests.complete_run.automation"
        if component == "controller"
        else "tests.complete_run.reporter"
    ) in marker.read_text(encoding="utf-8")
    if component == "controller":
        assert "--account\ne3sm\n" in marker.read_text(encoding="utf-8")
        assert "--qos\nregular\n" in marker.read_text(encoding="utf-8")
    assert (
        tmp_path / "results" / "automation" / "controller-environment.lock"
    ).exists()
    # The same manual wrapper still refuses to operate while an update holds the lock.
    marker.unlink()
    with scrontab._controller_lock(tmp_path / "results"):
        result = subprocess.run(
            [str(script), str(config), "--manual"],
            env=environment,
            check=False,
            capture_output=True,
        )
        assert result.returncode == 75
    assert not marker.exists()


def test_make_surface_and_environment_argument_forwarding(tmp_path: Path) -> None:
    makefile = (ops._CHECKOUT / "Makefile").read_text(encoding="utf-8")
    for prefix in (
        "complete-run-" + "ops-",
        "complete-run-" + "scron-",
        "ops-env-create",
        "ops-env-update",
        "ops-env-show",
        "ops-schedule-",
    ):
        assert prefix not in makefile
    result = subprocess.run(
        ["make", "--no-print-directory", "ops-help"],
        cwd=ops._CHECKOUT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Everyday commands" in result.stdout
    config = _config(tmp_path)
    result = subprocess.run(
        ["make", "--no-print-directory", "ops-logs", f"CONFIG={config}", "LINES=1"],
        cwd=ops._CHECKOUT,
        capture_output=True,
        text=True,
        check=True,
    )
    assert str(config) in result.stdout


def test_token_creation_never_prints_secret(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]
) -> None:
    config = _config(tmp_path)
    monkeypatch.setattr(ops.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(ops.getpass, "getpass", lambda prompt: "secret-value")
    ops._token_create(config, None)
    assert (tmp_path / "token").stat().st_mode & 0o777 == 0o600
    assert "secret-value" not in capsys.readouterr().out
    with pytest.raises(FileExistsError):
        ops._token_create(config, None)


def test_disable_needs_no_configuration(monkeypatch: pytest.MonkeyPatch) -> None:
    calls = []
    monkeypatch.setattr(scrontab, "remove_scrontab", lambda: calls.append("remove"))
    assert (
        ops.main(["disable", "--config", "/missing/controller.env", "--confirm", "YES"])
        == 0
    )
    assert calls == ["remove"]


def test_make_ignores_ambient_confirm_config_and_lines(tmp_path: Path) -> None:
    config = _config(tmp_path)
    environment = {
        **os.environ,
        "CONFIRM": "YES",
        "CONFIG": "/missing/controller.env",
        "LINES": "invalid",
        "E3SM_DIAGS_OPS_CONFIG": str(config),
    }
    # Confirmation must fail before any wrapper or mutation can run.
    result = subprocess.run(
        ["make", "ops-run"],
        cwd=ops._CHECKOUT,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    assert result.returncode != 0
    assert "CONFIRM=YES" in result.stderr
    result = subprocess.run(
        ["make", "ops-logs"],
        cwd=ops._CHECKOUT,
        env=environment,
        text=True,
        capture_output=True,
        check=True,
    )
    assert str(config) in result.stdout
    assert "invalid" not in result.stderr


def test_failed_update_reports_previous_revision_and_repair_steps(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    config = _config(tmp_path)
    monkeypatch.setattr(ops, "_fast_forward", lambda path: "previous-sha")

    def fail(path: Path) -> None:
        raise FileNotFoundError("controller environment missing")

    monkeypatch.setattr(ops, "validate_deployment", fail)
    with pytest.raises(
        RuntimeError, match="previous-sha.*deployment validation failed"
    ) as error:
        ops.update_checkout(config)
    assert "ops-disable CONFIRM=YES" in str(error.value)
    assert "checkout remains updated" in str(error.value)


def test_unfinished_git_operation_is_rejected(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    marker = tmp_path / "MERGE_HEAD"
    marker.touch()

    def git(arguments: list[str], **kwargs: object) -> subprocess.CompletedProcess[str]:
        command = arguments[3:]
        if command[0] == "status":
            output = ""
        elif "--git-path" in command:
            output = str(marker)
        else:
            output = "main"
        return subprocess.CompletedProcess(arguments, 0, output, "")

    monkeypatch.setattr(ops.subprocess, "run", git)
    with pytest.raises(ValueError, match="unfinished Git operation"):
        ops._fast_forward(tmp_path)


def test_token_with_explicit_path_needs_no_configuration(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    monkeypatch.setattr(ops.sys.stdin, "isatty", lambda: True)
    monkeypatch.setattr(ops.getpass, "getpass", lambda prompt: "secret-value")
    path = tmp_path / "new-token"
    ops.main(
        [
            "token-create",
            "--token-file",
            str(path),
            "--config",
            "/missing/controller.env",
        ]
    )
    assert path.stat().st_mode & 0o777 == 0o600
