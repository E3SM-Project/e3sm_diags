"""Tests for deferred complete-run reporting."""

from __future__ import annotations

import argparse
import json
import subprocess
from datetime import datetime, timedelta, timezone
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


def test_reporter_uses_sacct_when_squeue_reports_departed_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(args.results_root / "automation" / "run")

    def command(command: list[str]) -> str:
        if command[0] == "squeue":
            raise subprocess.CalledProcessError(
                1, command, stderr="slurm_load_jobs error: Invalid job id specified"
            )
        return "TIMEOUT|"

    monkeypatch.setattr(reporter, "_command", command)
    reporter.report_runs(args)

    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == "timed_out"


@pytest.mark.parametrize(
    "stage",
    [
        "submission_failed",
        "environment_failed",
        "diagnostics_failed",
        "cancelled",
        "timed_out",
        "slurm_failed",
        "job_completed_without_status",
    ],
)
def test_reporter_publishes_terminal_operational_failure(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, stage: str
):
    args = _args(tmp_path)
    run_root = args.results_root / "automation" / "run"
    _status(run_root, stage)
    publications: list[Path] = []

    def publish(_markdown: Path, receipt: Path, **_: object) -> dict[str, str]:
        publications.append(receipt)
        receipt.write_text('{"status": "published"}\n', encoding="utf-8")
        return {"status": "published"}

    monkeypatch.setattr(reporter, "publish_discussion", publish)

    reporter.report_runs(args)
    reporter.report_runs(args)

    assert publications == [run_root / "publication-receipt.json"]
    assert (
        "@E3SM-Project/e3sm-diags-admins: please review this operational failure"
        in (run_root / "automation-report.md").read_text(encoding="utf-8")
    )


def test_reporter_publishes_stalled_active_job(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    args.stall_threshold_hours = 24
    run_root = args.results_root / "automation" / "run"
    status = _status(run_root)
    payload = json.loads(status.read_text(encoding="utf-8"))
    payload["submitted_at_utc"] = (
        datetime.now(timezone.utc) - timedelta(hours=25)
    ).isoformat()
    status.write_text(json.dumps(payload), encoding="utf-8")
    monkeypatch.setattr(reporter, "_command", lambda _: "123 RUNNING")
    published: list[Path] = []

    def publish(_markdown: Path, receipt: Path, **_: object) -> dict[str, str]:
        published.append(receipt)
        receipt.write_text('{"status": "published"}\n', encoding="utf-8")
        return {"status": "published"}

    monkeypatch.setattr(reporter, "publish_discussion", publish)

    reporter.report_runs(args)

    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == "stalled"
    assert published == [run_root / "publication-receipt.json"]


def test_reporter_reuses_legacy_comparison_receipt(tmp_path: Path):
    run_root = tmp_path / "run"
    comparison = run_root / "comparison" / "set" / "comparison-report.json"
    comparison.parent.mkdir(parents=True)
    comparison.write_text("{}\n", encoding="utf-8")
    legacy_receipt = comparison.parent / "publication-receipt.json"
    legacy_receipt.write_text("{}\n", encoding="utf-8")

    assert reporter._publication_receipt_path(run_root) == legacy_receipt


def test_reporter_recovers_indeterminate_submission_handoff(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    status = _status(
        args.results_root / "automation" / "run", "submission_handoff_indeterminate"
    )
    monkeypatch.setattr(
        reporter,
        "_command",
        lambda command: "" if command[0] == "squeue" else "COMPLETED|",
    )

    reporter.report_runs(args)

    assert json.loads(status.read_text(encoding="utf-8"))["stage"] == (
        "job_completed_without_status"
    )


def test_reporter_does_not_repeat_final_report(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    args = _args(tmp_path)
    run_root = args.results_root / "automation" / "run"
    _status(run_root, "passed")
    reporter.report_runs(args)
    monkeypatch.setattr(reporter, "render_report", pytest.fail)

    reporter.report_runs(args)


def test_reporter_uses_report_title_for_publication(tmp_path: Path):
    run_root = tmp_path / "run"
    run_root.mkdir()
    (run_root / "automation-report.json").write_text(
        json.dumps(
            {"title": "E3SM Diags complete-run report — abc — 2026-09-25 17:39 UTC"}
        ),
        encoding="utf-8",
    )

    assert reporter._report_title(run_root) == (
        "E3SM Diags complete-run report — abc — 2026-09-25 17:39 UTC"
    )


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
    assert "controller-environment.lock" in controller
    assert "controller-environment.lock" in reporter_wrapper
    assert "TZ=America/Los_Angeles" in controller
    assert "TZ=America/Los_Angeles" in reporter_wrapper
