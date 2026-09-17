"""Tests for NERSC complete-run orchestration helpers."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pytest

from tests.complete_run import automation


def test_resolve_main_sha_fetches_then_resolves(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    calls: list[list[str]] = []

    def command(args: list[str], *, cwd: Path | None = None) -> str:
        calls.append(args)
        return "a" * 40 if args[1] == "rev-parse" else ""

    monkeypatch.setattr(automation, "_command", command)

    assert automation.resolve_main_sha(tmp_path) == "a" * 40
    assert calls == [
        ["git", "fetch", "origin", "main"],
        ["git", "rev-parse", "origin/main"],
    ]
    assert (
        automation.environment_name("abcdef123456789") == "e3sm_diags_ci_abcdef123456"
    )


def test_job_script_runs_diagnostics_then_full_comparison(tmp_path: Path):
    script = automation._job_script(
        tmp_path / "worktree",
        tmp_path / "env",
        tmp_path / "result",
        tmp_path / "comparison",
        tmp_path / "status.json",
        "abc",
        ["lat_lon"],
    )

    assert "tests.complete_run.run" in script
    assert "--workflow-revision abc" in script
    assert "tests.complete_run.compare" in script
    assert "--write-diff-html" in script
    assert '"diagnostics_failed"' in script
    assert '"comparison_failed"' in script


@pytest.mark.parametrize(
    ("sacct_output", "expected"),
    [
        ("CANCELLED by 1|\n", "cancelled"),
        ("TIMEOUT|\n", "timed_out"),
        ("FAILED|\n", "slurm_failed"),
        ("COMPLETED|\n", "job_completed_without_status"),
    ],
)
def test_terminal_stage(
    monkeypatch: pytest.MonkeyPatch, sacct_output: str, expected: str
):
    monkeypatch.setattr(automation, "_command", lambda _: sacct_output)
    assert automation._terminal_stage("1") == expected


def test_submission_failure_writes_machine_readable_status(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    args = argparse.Namespace(
        repo=tmp_path,
        results_root=tmp_path / "results",
        worktree_root=tmp_path / "worktrees",
        environment_root=tmp_path / "envs",
        account="e3sm",
        qos="regular",
        walltime="01:00:00",
        poll_seconds=0,
        sets=["lat_lon"],
        cfs_root=tmp_path,
        portal_root="https://portal.example",
    )
    monkeypatch.setattr(automation, "resolve_main_sha", lambda _: "a" * 40)
    monkeypatch.setattr(
        automation,
        "_command",
        lambda *_, **__: (_ for _ in ()).throw(subprocess.CalledProcessError(1, "git")),
    )

    assert automation.run_automation(args) == 1
    statuses = list((tmp_path / "results" / "automation").glob("*/status.json"))
    assert len(statuses) == 1
    assert (
        json.loads(statuses[0].read_text(encoding="utf-8"))["stage"]
        == "submission_failed"
    )


def test_existing_sha_environment_is_reported_without_submitting(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    environment_root = tmp_path / "envs"
    (environment_root / automation.environment_name("a" * 40)).mkdir(parents=True)
    args = argparse.Namespace(
        repo=tmp_path,
        results_root=tmp_path / "results",
        worktree_root=tmp_path / "worktrees",
        environment_root=environment_root,
        account="e3sm",
        qos="regular",
        walltime="01:00:00",
        poll_seconds=0,
        sets=["lat_lon"],
        cfs_root=tmp_path,
        portal_root="https://portal.example",
    )
    monkeypatch.setattr(automation, "resolve_main_sha", lambda _: "a" * 40)
    monkeypatch.setattr(
        automation, "_command", lambda *_, **__: pytest.fail("command run")
    )

    assert automation.run_automation(args) == 1
    status = json.loads(
        next((tmp_path / "results" / "automation").glob("*/status.json")).read_text()
    )
    assert "already exists" in status["error"]
