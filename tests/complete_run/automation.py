"""Orchestrate an immutable, NERSC complete-run environment regression."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from tests.complete_run.params import DEFAULT_RESULTS_DIR
from tests.complete_run.report import (
    DEFAULT_CFS_ROOT,
    DEFAULT_PORTAL_ROOT,
    render_report,
    write_report,
)
from tests.complete_run.run import DEFAULT_SETS_TO_RUN


def _command(args: list[str], *, cwd: Path | None = None) -> str:
    """Run a command and return stdout, raising on an execution failure."""
    return subprocess.run(
        args, check=True, capture_output=True, text=True, cwd=cwd
    ).stdout.strip()


def resolve_main_sha(repo: Path) -> str:
    """Fetch and resolve origin/main to one immutable commit SHA."""
    _command(["git", "fetch", "origin", "main"], cwd=repo)
    return _command(["git", "rev-parse", "origin/main"], cwd=repo)


def environment_name(sha: str) -> str:
    """Return the SHA-qualified Conda environment name."""
    return f"e3sm_diags_ci_{sha[:12]}"


def _write_json(path: Path, payload: dict[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _job_script(
    worktree: Path,
    prefix: Path,
    result_dir: Path,
    comparison_dir: Path,
    status_path: Path,
    sha: str,
    selected_sets: list[str],
) -> str:
    run_command = [
        "conda",
        "run",
        "-p",
        str(prefix),
        "python",
        "-m",
        "tests.complete_run.run",
        "--results-dir",
        str(result_dir),
        "--workflow-revision",
        sha,
    ]
    for diagnostic_set in selected_sets:
        run_command.extend(["--set", diagnostic_set])
    compare_command = [
        "conda",
        "run",
        "-p",
        str(prefix),
        "python",
        "-m",
        "tests.complete_run.compare",
        "--dev-dir",
        str(result_dir),
        "--baseline-dir",
        str(Path(DEFAULT_RESULTS_DIR) / "latest-main"),
        "--report-dir",
        str(comparison_dir),
        "--write-diff-pngs",
        "--write-diff-html",
    ]
    status = shlex.quote(str(status_path))
    return "\n".join(
        (
            "#!/bin/bash",
            "set -uo pipefail",
            f"cd {shlex.quote(str(worktree))}",
            f"if ! {' '.join(map(shlex.quote, run_command))}; then",
            f"  printf '%s\\n' '{{\"stage\": \"diagnostics_failed\"}}' > {status}",
            "  exit 0",
            "fi",
            f"if ! {' '.join(map(shlex.quote, compare_command))}; then",
            f"  printf '%s\\n' '{{\"stage\": \"comparison_failed\"}}' > {status}",
            "  exit 0",
            "fi",
            f"printf '%s\\n' '{{\"stage\": \"passed\"}}' > {status}",
            "",
        )
    )


def _terminal_stage(job_id: str) -> str:
    state = (
        _command(["sacct", "-j", job_id, "--format=State", "--noheader", "--parsable2"])
        .splitlines()[0]
        .split("|")[0]
    )
    if state.startswith("CANCELLED"):
        return "cancelled"
    if state.startswith("TIMEOUT"):
        return "timed_out"
    return "slurm_failed" if not state.startswith("COMPLETED") else "incomplete"


def run_automation(args: argparse.Namespace) -> int:
    """Create the environment, submit the job, monitor it, and render a report."""
    repo = args.repo.resolve()
    sha = resolve_main_sha(repo)
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    run_root = args.results_root / "automation" / f"{sha[:12]}-{stamp}"
    worktree = args.worktree_root / f"complete-run-{sha[:12]}-{stamp}"
    prefix = args.environment_root / environment_name(sha)
    result_dir = args.results_root / f"main-{sha[:12]}-{stamp}"
    comparison_dir = run_root / "comparison"
    status_path = run_root / "status.json"
    selected_sets = args.sets or DEFAULT_SETS_TO_RUN
    status: dict[str, object] = {
        "stage": "submission_failed",
        "git_sha": sha,
        "selected_sets": selected_sets,
        "environment_name": environment_name(sha),
        "environment_prefix": str(prefix),
        "result_dir": str(result_dir),
    }
    try:
        _command(["git", "worktree", "add", "--detach", str(worktree), sha], cwd=repo)
        _command(
            [
                "conda",
                "env",
                "create",
                "--prefix",
                str(prefix),
                "--file",
                str(worktree / "conda-env" / "ci.yml"),
            ]
        )
        _command(
            ["conda", "run", "-p", str(prefix), "pip", "install", "."], cwd=worktree
        )
        script_path = run_root / "complete-run.sbatch"
        script_path.parent.mkdir(parents=True, exist_ok=True)
        script_path.write_text(
            _job_script(
                worktree,
                prefix,
                result_dir,
                comparison_dir,
                status_path,
                sha,
                selected_sets,
            ),
            encoding="utf-8",
        )
        job_id = _command(
            [
                "sbatch",
                "--parsable",
                "--account",
                args.account,
                "--qos",
                args.qos,
                "--constraint",
                "cpu",
                "--nodes",
                "1",
                "--time",
                args.walltime,
                "--output",
                str(run_root / "slurm-%j.out"),
                str(script_path),
            ]
        ).split(";", 1)[0]
        status["job_id"] = job_id
        _write_json(status_path, status)
        while _command(["squeue", "-h", "-j", job_id]):
            time.sleep(args.poll_seconds)
        job_status = _load_job_status(status_path, status)
        if job_status["stage"] == "submission_failed":
            job_status["stage"] = _terminal_stage(job_id)
        _write_json(status_path, job_status)
    except (OSError, subprocess.CalledProcessError, IndexError) as error:
        status["error"] = str(error)
        _write_json(status_path, status)
    comparison_report = (
        next(comparison_dir.glob("*/comparison-report.json"), None)
        if comparison_dir.is_dir()
        else None
    )
    report = render_report(
        status_path,
        comparison_report,
        cfs_root=args.cfs_root,
        portal_root=args.portal_root,
    )
    write_report(report, run_root)
    return 0 if report["status"] == "passed" else 1


def _load_job_status(path: Path, fallback: dict[str, object]) -> dict[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return fallback
    return {**fallback, **payload}


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--repo", type=Path, default=Path.cwd())
    parser.add_argument("--results-root", type=Path, default=Path(DEFAULT_RESULTS_DIR))
    parser.add_argument("--worktree-root", type=Path, required=True)
    parser.add_argument("--environment-root", type=Path, required=True)
    parser.add_argument("--account", required=True)
    parser.add_argument("--qos", default="regular")
    parser.add_argument("--walltime", default="01:00:00")
    parser.add_argument("--poll-seconds", type=int, default=60)
    parser.add_argument(
        "--set", dest="sets", action="append", choices=DEFAULT_SETS_TO_RUN
    )
    parser.add_argument("--cfs-root", type=Path, default=DEFAULT_CFS_ROOT)
    parser.add_argument("--portal-root", default=DEFAULT_PORTAL_ROOT)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    """Run the NERSC login-node orchestration CLI."""
    return run_automation(_build_parser().parse_args(argv))


if __name__ == "__main__":
    raise SystemExit(main())
