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
        ["git", "fetch", "origin", "main:refs/remotes/origin/main"],
        ["git", "rev-parse", "origin/main"],
    ]
    assert (
        automation.environment_name("abcdef123456789", "20260917-120000")
        == "e3sm_diags_ci_abcdef123456_20260917-120000"
    )


def test_job_script_runs_diagnostics_then_full_comparison(tmp_path: Path):
    paths = {
        "worktree": tmp_path / "worktree",
        "prefix": tmp_path / "env",
        "result": tmp_path / "result",
        "comparison": tmp_path / "comparison",
        "status": tmp_path / "status.json",
    }
    script = automation._job_script(paths, "abc", ["lat_lon"])

    assert "tests.complete_run.run" in script
    assert "--workflow-revision abc" in script
    assert "tests.complete_run.compare" in script
    assert "--write-diff-html" in script
    assert "conda env create" in script
    assert "pip install ." in script
    assert '"environment_failed"' in script
    assert '"diagnostics_failed"' in script
    assert '"comparison_failed"' in script
    assert f"rm -rf {paths['prefix']}" in script
    assert script.index("conda env create") < script.index("pip install .")
    assert script.index("pip install .") < script.index("tests.complete_run.run")
    assert script.index("tests.complete_run.run") < script.index(
        "tests.complete_run.compare"
    )


def test_prepare_worktree_does_not_create_environment(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    calls: list[list[str]] = []

    def command(args: list[str], **_: object) -> str:
        calls.append(args)
        return ""

    monkeypatch.setattr(automation, "_command", command)
    paths = {"worktree": tmp_path / "worktree"}

    automation._prepare_worktree(tmp_path, paths, "a" * 40)

    assert calls == [
        ["git", "worktree", "add", "--detach", str(paths["worktree"]), "a" * 40]
    ]


def test_submit_job_uses_configured_slurm_resources(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    calls: list[list[str]] = []

    def command(args: list[str]) -> str:
        calls.append(args)
        return "123;cluster"

    monkeypatch.setattr(automation, "_command", command)
    args = argparse.Namespace(
        account="e3sm",
        qos="regular",
        constraint="cpu",
        nodes=1,
        walltime="02:00:00",
    )

    assert automation._submit_job(args, tmp_path / "job.sbatch", tmp_path) == "123"
    assert calls == [
        [
            "sbatch",
            "--parsable",
            "--account",
            "e3sm",
            "--qos",
            "regular",
            "--constraint",
            "cpu",
            "--nodes",
            "1",
            "--time",
            "02:00:00",
            "--output",
            str(tmp_path / "slurm-%j.out"),
            str(tmp_path / "job.sbatch"),
        ]
    ]


def test_submitted_job_status_returns_without_polling(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    paths = {
        "run_root": tmp_path / "run",
        "worktree": tmp_path / "worktree",
        "prefix": tmp_path / "environment",
        "result": tmp_path / "result",
        "comparison": tmp_path / "comparison",
        "status": tmp_path / "run" / "status.json",
    }
    status = automation._initial_status("a" * 40, paths, ["lat_lon"])
    args = argparse.Namespace()
    monkeypatch.setattr(automation, "_submit_job", lambda *_: "123")

    automation._submit_job_for_run(args, paths, "a" * 40, ["lat_lon"], status)

    payload = json.loads(paths["status"].read_text(encoding="utf-8"))
    assert payload["stage"] == "submitted"
    assert "submitted_at_utc" in payload


def test_run_paths_keep_results_with_immutable_run_root(tmp_path: Path):
    args = argparse.Namespace(
        results_root=tmp_path / "results",
        worktree_root=tmp_path / "worktrees",
        environment_root=tmp_path / "environments",
    )

    paths = automation._build_run_paths(args, "a" * 40)

    assert paths["result"].parent == paths["run_root"]
    assert paths["result"].name == "results"


def test_parser_parses_configured_node_count(tmp_path: Path):
    args = automation._build_parser().parse_args(
        [
            "--worktree-root",
            str(tmp_path / "worktrees"),
            "--environment-root",
            str(tmp_path / "environments"),
            "--account",
            "e3sm",
            "--nodes",
            "2",
        ]
    )

    assert args.nodes == 2


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
        sets=["lat_lon"],
        cfs_root=tmp_path,
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


def test_post_submission_status_failure_preserves_worktree(
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
        constraint="cpu",
        sets=["lat_lon"],
    )
    monkeypatch.setattr(automation, "resolve_main_sha", lambda _: "a" * 40)
    monkeypatch.setattr(automation, "_prepare_worktree", lambda *_: None)
    monkeypatch.setattr(automation, "_submit_job", lambda *_: "123")
    original_write = automation._write_json
    writes = 0

    def write_json(path: Path, payload: dict[str, object]) -> None:
        nonlocal writes
        writes += 1
        if writes == 2:
            raise OSError("CFS unavailable")
        original_write(path, payload)

    removed: list[Path] = []
    monkeypatch.setattr(automation, "_write_json", write_json)
    monkeypatch.setattr(
        automation, "_remove_worktree", lambda _, path: removed.append(path)
    )

    assert automation.run_automation(args) == 1
    assert removed == []
    status = next((tmp_path / "results" / "automation").glob("*/status.json"))
    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == (
        "submission_handoff_indeterminate"
    )


def test_revision_resolution_failure_writes_completion_metadata(
    monkeypatch: pytest.MonkeyPatch, tmp_path: Path
):
    completion_file = tmp_path / "completion.json"
    args = argparse.Namespace(
        repo=tmp_path,
        results_root=tmp_path / "results",
        worktree_root=tmp_path / "worktrees",
        environment_root=tmp_path / "envs",
        account="e3sm",
        qos="regular",
        walltime="01:00:00",
        sets=["lat_lon"],
        cfs_root=tmp_path,
        completion_file=completion_file,
    )
    error = subprocess.CalledProcessError(
        128, ["git", "fetch", "origin", "main"], stderr="fatal: authentication failed"
    )
    monkeypatch.setattr(
        automation,
        "resolve_main_sha",
        lambda _: (_ for _ in ()).throw(error),
    )

    assert automation.run_automation(args) == 1
    run_root = Path(json.loads(completion_file.read_text(encoding="utf-8"))["run_root"])
    status = json.loads((run_root / "status.json").read_text(encoding="utf-8"))
    assert status["stage"] == "submission_failed"
    assert status["git_sha"] is None
    assert "fatal: authentication failed" in status["error"]
    assert not (run_root / "automation-report.json").exists()
