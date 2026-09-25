"""Orchestrate an immutable, NERSC complete-run environment regression."""

from __future__ import annotations

import argparse
import json
import shlex
import subprocess
from datetime import datetime, timezone
from pathlib import Path
from typing import Sequence

from e3sm_diags.logger import _setup_child_logger
from tests.complete_run.params import DEFAULT_RESULTS_DIR
from tests.complete_run.run import DEFAULT_SETS_TO_RUN

logger = _setup_child_logger(__name__)


def main(argv: Sequence[str] | None = None) -> int:
    """Run the NERSC login-node orchestration CLI."""
    args = _build_parser().parse_args(argv)

    return run_automation(args)


def resolve_main_sha(repo: Path) -> str:
    """Fetch and resolve ``origin/main`` to one immutable commit SHA."""
    # An explicitly named fetch updates FETCH_HEAD, but it does not necessarily
    # create origin/main (for example, in a single-branch controller checkout).
    # Map the remote branch explicitly so the revision resolved below always
    # exists in the local repository.
    _command(["git", "fetch", "origin", "main:refs/remotes/origin/main"], cwd=repo)
    return _command(["git", "rev-parse", "origin/main"], cwd=repo)


def environment_name(sha: str, run_id: str) -> str:
    """Return a fresh Conda environment name qualified by SHA and run ID."""
    return f"e3sm_diags_ci_{sha[:12]}_{run_id}"


def run_automation(args: argparse.Namespace) -> int:
    """Create a worktree and submit diagnostics without waiting for Slurm."""
    repo = args.repo.resolve()
    selected_sets = args.sets or DEFAULT_SETS_TO_RUN
    try:
        sha = resolve_main_sha(repo)
    except (OSError, subprocess.CalledProcessError) as error:
        paths = _build_run_paths(args, "unresolved")
        status = _initial_status("unresolved", paths, selected_sets)
        status["stage"] = "submission_failed"
        status["git_sha"] = None
        status["error"] = _command_error(error)
        _write_json(paths["status"], status)
        _write_completion_file(args, paths["run_root"])
        return 1

    paths = _build_run_paths(args, sha)
    status = _initial_status(sha, paths, selected_sets)
    _write_json(paths["status"], status)

    try:
        _prepare_worktree(repo, paths, sha)
        _submit_job_for_run(args, paths, sha, selected_sets, status)
    except (OSError, subprocess.CalledProcessError, IndexError) as error:
        status["error"] = _command_error(error)
        # A job ID means Slurm may already be using this worktree.  In
        # particular, CFS can fail between ``sbatch`` returning successfully
        # and the submitted status being persisted.  Never invalidate that
        # allocation by removing its checkout.
        if "job_id" in status:
            status["stage"] = "submission_handoff_indeterminate"
            logger.warning(
                "Submitted complete-run job %s could not be fully recorded: %s",
                status["job_id"],
                error,
            )
        else:
            status["stage"] = "submission_failed"
            _remove_worktree(repo, paths["worktree"])
        try:
            _write_json(paths["status"], status)
        except OSError as write_error:
            logger.warning(
                "Unable to persist complete-run submission recovery status in %s: %s",
                paths["run_root"],
                write_error,
            )

    _write_completion_file(args, paths["run_root"])
    if status.get("stage") == "submitted":
        return 0

    return 1


def _build_run_paths(args: argparse.Namespace, sha: str) -> dict[str, Path]:
    """Build immutable paths for one orchestration attempt."""
    stamp = datetime.now(timezone.utc).strftime("%Y%m%d-%H%M%S")
    run_root = args.results_root / "automation" / f"{sha[:12]}-{stamp}"
    return {
        "run_root": run_root,
        "worktree": args.worktree_root / f"complete-run-{sha[:12]}-{stamp}",
        "prefix": args.environment_root / environment_name(sha, stamp),
        "result": run_root / "results",
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


def _prepare_worktree(repo: Path, paths: dict[str, Path], sha: str) -> None:
    """Create the detached worktree before submitting the compute allocation."""
    worktree = paths["worktree"]
    _command(["git", "worktree", "add", "--detach", str(worktree), sha], cwd=repo)


def _remove_worktree(repo: Path, worktree: Path) -> None:
    """Best-effort cleanup after a submission failure before Slurm owns the run."""
    if worktree.exists():
        try:
            _command(["git", "worktree", "remove", "--force", str(worktree)], cwd=repo)
        except (OSError, subprocess.CalledProcessError) as error:
            logger.warning(
                "Unable to remove failed complete-run worktree %s: %s", worktree, error
            )


def _submit_job_for_run(
    args: argparse.Namespace,
    paths: dict[str, Path],
    sha: str,
    selected_sets: list[str],
    status: dict[str, object],
) -> None:
    """Submit the compute job and preserve the handoff status for a reporter."""
    script_path = paths["run_root"] / "complete-run.sbatch"
    script_path.parent.mkdir(parents=True, exist_ok=True)
    script_path.write_text(_job_script(paths, sha, selected_sets), encoding="utf-8")
    job_id = _submit_job(args, script_path, paths["run_root"])
    status["job_id"] = job_id
    status["stage"] = "submitted"
    _write_json(paths["status"], status)


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
            args.constraint,
            "--nodes",
            str(args.nodes),
            "--time",
            args.walltime,
            "--output",
            str(run_root / "slurm-%j.out"),
            str(script_path),
        ]
    )
    return output.split(";", 1)[0]


def _write_completion_file(args: argparse.Namespace, run_root: Path) -> None:
    """Persist the run root for a controller that invokes this CLI."""
    completion_file = getattr(args, "completion_file", None)
    if completion_file is not None:
        _write_json(completion_file, {"run_root": str(run_root)})


def _job_script(paths: dict[str, Path], sha: str, selected_sets: list[str]) -> str:
    """Build the batch script that provisions and runs the diagnostics environment."""
    create_environment_command = [
        "conda",
        "env",
        "create",
        "--prefix",
        str(paths["prefix"]),
        "--file",
        str(paths["worktree"] / "conda-env" / "ci.yml"),
    ]
    install_command = [
        "conda",
        "run",
        "-p",
        str(paths["prefix"]),
        "pip",
        "install",
        ".",
    ]
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
    prefix = shlex.quote(str(paths["prefix"]))

    def write_status(stage: str) -> str:
        payload = {
            "stage": stage,
            "git_sha": sha,
            "selected_sets": selected_sets,
            "environment_name": paths["prefix"].name,
            "environment_prefix": str(paths["prefix"]),
            "result_dir": str(paths["result"]),
        }
        return f"printf '%s\\n' {shlex.quote(json.dumps(payload))} > {status}"

    return "\n".join(
        (
            "#!/bin/bash",
            "set -euo pipefail",
            f"cd {shlex.quote(str(paths['worktree']))}",
            f"if ! {' '.join(map(shlex.quote, create_environment_command))}; then",
            f"  rm -rf {prefix}",
            f"  {write_status('environment_failed')}",
            "  exit 0",
            "fi",
            f"if ! {' '.join(map(shlex.quote, install_command))}; then",
            f"  rm -rf {prefix}",
            f"  {write_status('environment_failed')}",
            "  exit 0",
            "fi",
            f"if ! {' '.join(map(shlex.quote, run_command))}; then",
            f"  {write_status('diagnostics_failed')}",
            "  exit 0",
            "fi",
            f"if ! {' '.join(map(shlex.quote, compare_command))}; then",
            f"  {write_status('comparison_failed')}",
            "  exit 0",
            "fi",
            write_status("passed"),
            "",
        )
    )


def _command(args: list[str], *, cwd: Path | None = None) -> str:
    """Run a command and return stripped standard output."""
    completed = subprocess.run(
        args, check=True, capture_output=True, text=True, cwd=cwd
    )
    return completed.stdout.strip()


def _command_error(error: OSError | subprocess.CalledProcessError | IndexError) -> str:
    """Return an actionable error message, including failed-command stderr."""
    if isinstance(error, subprocess.CalledProcessError):
        stderr = error.stderr.strip() if isinstance(error.stderr, str) else ""
        if stderr:
            return f"{error}: {stderr}"

    return str(error)


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
    parser.add_argument("--nodes", type=int, default=1)
    parser.add_argument("--walltime", default="02:00:00")
    parser.add_argument("--constraint", default="cpu")
    parser.add_argument(
        "--set", dest="sets", action="append", choices=DEFAULT_SETS_TO_RUN
    )
    parser.add_argument(
        "--completion-file",
        type=Path,
        help="Machine-readable location of the finished orchestration artifacts.",
    )
    return parser


if __name__ == "__main__":
    raise SystemExit(main())
