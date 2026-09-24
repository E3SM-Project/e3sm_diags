"""Tests for deferred complete-run reporting."""

from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

import pytest

from tests.complete_run import reporter

_COMPLETE_RUN_ROOT = Path(__file__).parents[1] / "complete_run"


def _args(tmp_path: Path) -> argparse.Namespace:
    return argparse.Namespace(
        results_root=tmp_path / "results",
        repository_id="repository",
        category_id="category",
        token_file=tmp_path / "token",
        cfs_root=tmp_path,
        portal_root="https://portal.example",
    )


def _status(run_root: Path, stage: str = "submitted") -> Path:
    path = run_root / "status.json"
    path.parent.mkdir(parents=True)
    path.write_text(
        json.dumps(
            {"stage": stage, "job_id": "123", "result_dir": str(run_root / "result")}
        ),
        encoding="utf-8",
    )
    return path


def test_reporter_leaves_queued_job_unchanged(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(args.results_root / "automation" / "run")
    monkeypatch.setattr(reporter, "_command", lambda _: "123 queued")

    assert reporter.report_runs(args) == 0
    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == "submitted"


def test_reporter_renders_batch_written_status(tmp_path: Path):
    args = _args(tmp_path)
    run_root = args.results_root / "automation" / "run"
    _status(run_root, "passed")

    reporter.report_runs(args)

    assert (run_root / "automation-report.json").is_file()


def test_reporter_classifies_departed_job_without_status(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(args.results_root / "automation" / "run")
    monkeypatch.setattr(
        reporter,
        "_command",
        lambda command: "" if command[0] == "squeue" else "COMPLETED|",
    )

    reporter.report_runs(args)

    assert (
        json.loads(status.read_text(encoding="utf-8"))["stage"]
        == "job_completed_without_status"
    )


def test_reporter_retries_when_accounting_lags(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(args.results_root / "automation" / "run")

    def command(command: list[str]) -> str:
        if command[0] == "squeue":
            return ""
        raise subprocess.CalledProcessError(1, command)

    monkeypatch.setattr(reporter, "_command", command)
    reporter.report_runs(args)

    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == "submitted"
    assert not (status.parent / "automation-report.json").exists()


def test_reporter_leaves_run_for_a_later_tick_when_squeue_fails(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(args.results_root / "automation" / "run")
    monkeypatch.setattr(
        reporter,
        "_command",
        lambda command: (_ for _ in ()).throw(
            subprocess.CalledProcessError(1, command)
        ),
    )

    reporter.report_runs(args)

    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == "submitted"


def test_reporter_does_not_repeat_final_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    run_root = args.results_root / "automation" / "run"
    _status(run_root, "passed")
    reporter.report_runs(args)
    monkeypatch.setattr(reporter, "render_report", pytest.fail)

    reporter.report_runs(args)


def test_shell_wrappers_invoke_their_separate_entry_points():
    controller = (_COMPLETE_RUN_ROOT / "complete-run-controller.sh").read_text(
        encoding="utf-8"
    )
    reporter_wrapper = (_COMPLETE_RUN_ROOT / "complete-run-reporter.sh").read_text(
        encoding="utf-8"
    )

    assert "tests.complete_run.automation" in controller
    assert "tests.complete_run.reporter" not in controller
    assert "tests.complete_run.reporter" in reporter_wrapper
    assert "reporter.lock" in reporter_wrapper
