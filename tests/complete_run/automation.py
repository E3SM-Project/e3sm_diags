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


def resolve_main_sha(repo: Path) -> str:
    """Fetch and resolve ``origin/main`` to one immutable commit SHA."""
    _command(["git", "fetch", "origin", "main"], cwd=repo)
    return _command(["git", "rev-parse", "origin/main"], cwd=repo)


def environment_name(sha: str, run_id: str) -> str:
    """Return a fresh Conda environment name qualified by SHA and run ID."""
    return f"e3sm_diags_ci_{sha[:12]}_{run_id}"


def run_automation(args: argparse.Namespace) -> int:
    """Create an environment, submit diagnostics, and render its report."""
    repo = args.repo.resolve()
    sha = resolve_main_sha(repo)
    paths = _build_run_paths(args, sha)
    selected_sets = args.sets or DEFAULT_SETS_TO_RUN
    status = _initial_status(sha, paths, selected_sets)

    try:
        _prepare_environment(repo, paths, sha)
        _submit_and_monitor_job(args, paths, sha, selected_sets, status)
    except (OSError, subprocess.CalledProcessError, IndexError) as error:
        status["error"] = str(error)
        _write_json(paths["status"], status)

    comparison_report = _comparison_report(paths["comparison"])
    report = render_report(
        paths["status"],
        comparison_report,
        cfs_root=args.cfs_root,
        portal_root=args.portal_root,
    )
    write_report(report, paths["run_root"])
    _write_completion_file(args, paths["run_root"])

    if report["status"] == "passed":
        return 0

    return 1


def main(argv: Sequence[str] | None = None) -> int:
    """Run the NERSC login-node orchestration CLI."""
    args = _build_parser().parse_args(argv)
    return run_automation(args)


def _build_run_paths(args: argparse.Namespace, sha: str) -> dict[str, Path]:
    """Build immutable paths for one orchestration attempt."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    run_root = args.results_root / "automation" / f"{sha[:12]}-{stamp}"
    return {
        "run_root": run_root,
        "worktree": args.worktree_root / f"complete-run-{sha[:12]}-{stamp}",
        "prefix": args.environment_root / environment_name(sha, stamp),
        "result": args.results_root / f"main-{sha[:12]}-{stamp}",
        "comparison": run_root / "comparison",
        "status": run_root / "status.json",
    }


def _initial_status(
    sha: str, paths: dict[str, Path], selected_sets: list[str]
) -> dict[str, object]:
    """Build the status preserved for an unsuccessful submission."""
    return {
        "stage": "submission_failed",
        "git_sha": sha,
        "selected_sets": selected_sets,
        "environment_name": paths["prefix"].name,
        "environment_prefix": str(paths["prefix"]),
        "result_dir": str(paths["result"]),
    }


def _prepare_environment(repo: Path, paths: dict[str, Path], sha: str) -> None:
    """Create the detached worktree and fresh SHA-qualified environment."""
    worktree = paths["worktree"]
    _command(["git", "worktree", "add", "--detach", str(worktree), sha], cwd=repo)
    _command(
        [
            "conda",
            "env",
            "create",
            "--prefix",
            str(paths["prefix"]),
            "--file",
            str(worktree / "conda-env" / "ci.yml"),
        ]
    )
    _command(
        ["conda", "run", "-p", str(paths["prefix"]), "pip", "install", "."],
        cwd=worktree,
    )


def _submit_and_monitor_job(
    args: argparse.Namespace,
    paths: dict[str, Path],
    sha: str,
    selected_sets: list[str],
    status: dict[str, object],
) -> None:
    """Submit the compute job and preserve its final Slurm status."""
    script_path = paths["run_root"] / "complete-run.sbatch"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(_job_script(paths, sha, selected_sets), encoding="utf-8")
    job_id = _submit_job(args, script_path, paths["run_root"])
    status["job_id"] = job_id
    _write_json(paths["status"], status)

    while _command(["squeue", "-h", "-j", job_id]):
        time.sleep(args.poll_seconds)

    final_status = _load_job_status(paths["status"], status)
    if final_status["stage"] == "submission_failed":
        final_status["stage"] = _terminal_stage(job_id)
    _write_json(paths["status"], final_status)


def _submit_job(args: argparse.Namespace, script_path: Path, run_root: Path) -> str:
    """Submit a CPU Slurm job and return its parsable job ID."""
    output = _command(
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
    )
    return output.split(";", 1)[0]


def _comparison_report(comparison_dir: Path) -> Path | None:
    """Return the comparison report produced by this orchestration attempt."""
    if not comparison_dir.is_dir():
        return None

    return next(comparison_dir.glob("*/comparison-report.json"), None)


def _write_completion_file(args: argparse.Namespace, run_root: Path) -> None:
    """Persist the run root for a controller that invokes this CLI."""
    completion_file = getattr(args, "completion_file", None)
    if completion_file is not None:
        _write_json(completion_file, {"run_root": str(run_root)})


def _job_script(paths: dict[str, Path], sha: str, selected_sets: list[str]) -> str:
    """Build the batch script that records diagnostics and comparison outcomes."""
    run_command = [
        "conda",
        "run",
        "-p",
        str(paths["prefix"]),
        "python",
        "-m",
        "tests.complete_run.run",
        "--results-dir",
        str(paths["result"]),
        "--workflow-revision",
        sha,
    ]
    for diagnostic_set in selected_sets:
        run_command.extend(["--set", diagnostic_set])

    compare_command = [
        "conda",
        "run",
        "-p",
        str(paths["prefix"]),
        "python",
        "-m",
        "tests.complete_run.compare",
        "--dev-dir",
        str(paths["result"]),
        "--baseline-dir",
        str(Path(DEFAULT_RESULTS_DIR) / "latest-main"),
        "--report-dir",
        str(paths["comparison"]),
        "--write-diff-pngs",
        "--write-diff-html",
    ]
    status = shlex.quote(str(paths["status"]))
    return "\n".join(
        (
            "#!/bin/bash",
            "set -euo pipefail",
            f"cd {shlex.quote(str(paths['worktree']))}",
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
    """Classify a terminal Slurm state when the batch script wrote no status."""
    state = (
        _command(["sacct", "-j", job_id, "--format=State", "--noheader", "--parsable2"])
        .splitlines()[0]
        .split("|")[0]
    )
    if state.startswith("CANCELLED"):
        return "cancelled"
    if state.startswith("TIMEOUT"):
        return "timed_out"
    if state.startswith("COMPLETED"):
        return "job_completed_without_status"

    return "slurm_failed"


def _command(args: list[str], *, cwd: Path | None = None) -> str:
    """Run a command and return stripped standard output."""
    completed = subprocess.run(
        args, check=True, capture_output=True, text=True, cwd=cwd
    )
    return completed.stdout.strip()


def _load_job_status(path: Path, fallback: dict[str, object]) -> dict[str, object]:
    """Load an in-job status without discarding orchestration provenance."""
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return fallback

    return {**fallback, **payload}


def _write_json(path: Path, payload: dict[str, object]) -> None:
    """Write a machine-readable status or controller handoff file."""
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def _build_parser() -> argparse.ArgumentParser:
    """Build the explicitly invoked NERSC orchestration CLI."""
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
    parser.add_argument(
        "--completion-file",
        type=Path,
        help="Machine-readable location of the finished orchestration artifacts.",
    )
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
