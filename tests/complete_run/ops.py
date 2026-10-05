"""Operator interface for the persistent complete-run deployment."""

from __future__ import annotations

import argparse
import getpass
import json
import os
import shlex
import subprocess
import sys
from collections import deque
from pathlib import Path
from typing import Sequence

from tests.complete_run import scrontab

_CHECKOUT = Path(__file__).resolve().parents[2]
_EVERYDAY = ("status", "logs", "update", "run", "report", "help")
_ADMIN = ("init", "env", "enable", "disable", "token-create", "shortcut")
_HELP = """Everyday commands:
  make ops                         Read-only deployment dashboard
  make ops-logs [JOB=<id>] [LINES=100]
  make ops-update                  Fast-forward checkout and validate; no Conda update
  make ops-run CONFIRM=YES          Submit a catch-up run (no cron time/week guards)
  make ops-report CONFIRM=YES       Process completed runs and retry publication
  make ops-help

Administration:
  make ops-init OPERATIONS_DIR=/path [BRANCH=main] [REPOSITORY_URL=<url>]
  make ops-env ACTION=create
  make ops-env ACTION=update CONFIRM=YES
  make ops-enable CONFIRM=YES
  make ops-disable CONFIRM=YES
  make ops-token-create [TOKEN_FILE=/path]
  make ops-shortcut                 Print a shell function; do not install it

Configuration: CONFIG=/path/controller.env, then E3SM_DIAGS_OPS_CONFIG,
then controller.env in this checkout's parent. Setup never enables scheduling.
"""


def _output(text: str) -> None:
    """Write operator-facing output without configuring application logging."""
    sys.stdout.write(text + "\n")


def resolve_config(explicit: Path | None = None) -> Path:
    """Resolve configuration without searching the filesystem.

    Raises
    ------
    FileNotFoundError
        If the selected configuration does not exist.
    """
    selected = explicit or os.environ.get("E3SM_DIAGS_OPS_CONFIG")
    path = Path(selected) if selected else _CHECKOUT.parent / "controller.env"
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Controller configuration not found: {path}. Set CONFIG=/path/controller.env "
            "or E3SM_DIAGS_OPS_CONFIG; for initial setup use "
            "make ops-init OPERATIONS_DIR=/path and edit controller.env."
        )
    return path


def _inspect_command(arguments: list[str]) -> str:
    """Run a read-only tool, returning an actionable diagnostic on failure."""
    try:
        result = subprocess.run(arguments, capture_output=True, text=True, check=True)
        return result.stdout.strip() or "(none)"
    except FileNotFoundError:
        return f"Unavailable: {arguments[0]} is not installed or not on PATH."
    except subprocess.CalledProcessError as error:
        return f"Unavailable: {arguments[0]}: {(error.stderr or str(error)).strip()}"


def _metadata(path: Path) -> dict[str, object]:
    """Load bounded run metadata, retaining malformed-file diagnostics."""
    if not path.is_file():
        return {"unavailable": f"Not present: {path.name}"}
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
        if not isinstance(payload, dict):
            raise ValueError("expected a JSON object")
        return payload
    except (OSError, ValueError) as error:
        return {"unavailable": f"Malformed or unreadable metadata {path}: {error}"}


def dashboard(config_path: Path) -> None:
    """Show deployment health, schedule, jobs, and latest recorded run read-only."""
    config = scrontab._read_config(config_path)
    _output(f"Configuration: {config_path}")
    for key in (
        "REPOSITORY",
        "LOG_DIR",
        "CONDA_BASE",
        "CONTROLLER_ENV_PREFIX",
        "RESULTS_ROOT",
    ):
        value = config.get(key, "")
        _output(
            f"{key}: {value or '(missing)'} [present={bool(value) and Path(value).exists()}]"
        )
    try:
        scrontab.validate_config(config_path)
        _output("Configuration health: valid")
    except (OSError, ValueError) as error:
        _output(f"Configuration health: {error}")
    repository = config.get("REPOSITORY")
    if repository:
        for component in ("controller", "reporter"):
            script = (
                Path(repository)
                / "tests"
                / "complete_run"
                / f"complete-run-{component}.sh"
            )
            _output(
                f"Script: {script} [executable={script.is_file() and os.access(script, os.X_OK)}]"
            )
    _output("Installed schedule (future recurring occurrences, not previous outcomes):")
    try:
        _output(scrontab._read_scrontab() or "(no schedule)")
    except (OSError, subprocess.CalledProcessError) as error:
        _output(f"Schedule unavailable: {error}")
    _output("Slurm cron jobs / next eligible occurrences:")
    _output(
        _inspect_command(
            ["squeue", "--me", "-q", "cron", "-O", "JobID,Name,State,EligibleTime"]
        )
    )
    try:
        _latest_run(config)
    except OSError as error:
        _output(f"Latest automated run unavailable: {error}")


def _latest_run(config: dict[str, str]) -> None:
    """Inspect only immediate automated run directories, never results trees."""
    root_value = config.get("RESULTS_ROOT")
    if not root_value:
        _output("Latest automated run: RESULTS_ROOT is missing.")
        return
    root = Path(root_value) / "automation"
    candidates = sorted(
        (path for path in root.glob("*-????????-??????") if path.is_dir()),
        key=lambda path: path.name.rsplit("-", 2)[-2:],
    )
    if not candidates:
        _output(f"Latest automated run: none in {root}")
        return
    run = candidates[-1]
    _output(f"Latest automated run: {run}")
    status = _metadata(run / "status.json")
    _output("Recorded outcome: " + json.dumps(status, sort_keys=True))
    job = str(status.get("job_id", ""))
    if job.isdigit():
        _output("Recorded run Slurm accounting (not the next recurring occurrence):")
        _output(
            _inspect_command(
                ["sacct", "-j", job, "--format=JobID,State,ExitCode", "--noheader"]
            )
        )
    for filename in (
        "automation-report.json",
        "publication-receipt.json",
        "publication-failure.json",
    ):
        payload = _metadata(run / filename)
        if filename == "automation-report.json" and "unavailable" not in payload:
            payload = {key: payload.get(key) for key in ("status", "publication")}
        _output(f"{filename}: " + json.dumps(payload, sort_keys=True))
    _output(
        "Publication policy: qualifying failures only; a successful run need not have a Discussion."
    )


def logs(config_path: Path, *, job: str | None, lines: int) -> None:
    """Show bounded tails of recent controller and reporter logs.

    Raises
    ------
    ValueError
        If the job ID or line limit is invalid.
    """
    if lines < 1 or (job is not None and not job.isdigit()):
        raise ValueError(
            "LINES must be positive and JOB must be a numeric Slurm job ID."
        )
    config = scrontab._read_config(config_path)
    if not config.get("LOG_DIR"):
        raise ValueError("LOG_DIR is missing from controller configuration.")
    root = Path(config["LOG_DIR"])
    _output(f"Configuration: {config_path}\nLogs: {root}")
    for component in ("controller", "reporter"):
        paths = list(root.glob(f"complete-run-{component}-{job or '*'}.out"))
        paths = sorted(
            (path for path in paths if path.is_file()),
            key=lambda path: path.stat().st_mtime,
        )
        if not paths:
            _output(f"No {component} logs found" + (f" for JOB={job}." if job else "."))
        for path in paths[-1:]:
            _output(f"--- {path} (last {lines} lines) ---")
            with path.open(encoding="utf-8", errors="replace") as stream:
                sys.stdout.write("".join(deque(stream, maxlen=lines)))


def validate_deployment(config_path: Path) -> None:
    """Validate configuration, executable wrappers and controller imports.

    Raises
    ------
    FileNotFoundError
        If a required deployment artifact is absent or not executable.
    """
    scrontab.validate_config(config_path)
    environment = scrontab._controller_environment(config_path)
    for component in ("controller", "reporter"):
        script = (
            environment.repository
            / "tests"
            / "complete_run"
            / f"complete-run-{component}.sh"
        )
        if not script.is_file() or not os.access(script, os.X_OK):
            raise FileNotFoundError(
                f"Required executable script is missing or not executable: {script}"
            )
    if not environment.prefix.is_dir():
        raise FileNotFoundError(
            f"Controller environment is missing: {environment.prefix}; run make ops-env ACTION=create."
        )
    scrontab._verify_controller_environment(environment)


def update_checkout(config_path: Path) -> None:
    """Fast-forward a clean, attached checkout under the shared environment lock.

    Raises
    ------
    ValueError
        If the checkout is dirty, detached, untracked, or ahead/diverged.
    """
    environment = scrontab._controller_environment(config_path)
    with scrontab._controller_lock(environment.results_root):
        previous = _fast_forward(environment.repository)
        try:
            validate_deployment(config_path)
        except (OSError, ValueError, subprocess.CalledProcessError) as error:
            raise RuntimeError(
                f"Checkout updated from {previous} but deployment validation failed: {error}. "
                "The checkout remains updated. Disable scheduling with make ops-disable CONFIRM=YES "
                "before repair; inspect the checkout revision and, if dependencies changed, deliberately "
                "run make ops-env ACTION=update CONFIRM=YES, then ops-enable CONFIRM=YES. "
                f"Previous revision for manual recovery: {previous}."
            ) from error


def _fast_forward(repository: Path) -> str:
    """Reject unsafe Git states and merge only a fetched upstream fast-forward."""

    def git(*arguments: str) -> str:
        return subprocess.run(
            ["git", "-C", str(repository), *arguments],
            text=True,
            capture_output=True,
            check=True,
        ).stdout.strip()

    if git("status", "--porcelain"):
        raise ValueError(
            "Operations checkout is dirty; commit or remove local changes before ops-update."
        )
    branch = git("branch", "--show-current")
    if not branch:
        raise ValueError(
            "Operations checkout is detached; check out its deployment branch."
        )
    for marker in (
        "MERGE_HEAD",
        "CHERRY_PICK_HEAD",
        "REVERT_HEAD",
        "rebase-merge",
        "rebase-apply",
        "BISECT_LOG",
    ):
        if Path(
            git("rev-parse", "--path-format=absolute", "--git-path", marker)
        ).exists():
            raise ValueError(
                f"Operations checkout has an unfinished Git operation: {marker}"
            )
    upstream = git("rev-parse", "--abbrev-ref", "--symbolic-full-name", "@{upstream}")
    remote = git("config", "--get", f"branch.{branch}.remote")
    previous = git("rev-parse", "HEAD")
    git("fetch", remote)
    ahead, _ = git("rev-list", "--left-right", "--count", f"HEAD...{upstream}").split()
    if int(ahead):
        raise ValueError(
            "Operations checkout is ahead of or diverged from upstream; refusing update."
        )
    git("merge", "--ff-only", upstream)
    _output(
        f"Checkout: {repository}\nRevision: {previous} -> {git('rev-parse', 'HEAD')}"
    )
    return previous


def shortcut(config_path: Path) -> str:
    """Return a safely quoted Bash function usable outside the checkout."""
    config = scrontab._read_config(config_path)
    repository = config.get("REPOSITORY", "")
    if not repository or not Path(repository).is_absolute():
        raise ValueError("REPOSITORY must be an absolute configured checkout path.")
    return (
        "e3sm-ops() {\n"
        "  local target=ops\n"
        '  if [[ $# -gt 0 && "$1" != *=* ]]; then\n'
        '    case "$1" in\n'
        "      status) target=ops ;;\n"
        "      --help) target=ops-help ;;\n"
        '      logs|update|run|report|help|init|env|enable|disable|token-create|shortcut) target=ops-"$1" ;;\n'
        '      *) printf "%s\\n" "Unknown operations command: $1" >&2; return 2 ;;\n'
        "    esac\n    shift\n  fi\n"
        f'  command make -C {shlex.quote(repository)} "$target" "$@" CONFIG={shlex.quote(str(config_path))}\n'
        "}"
    )


def _token_create(config_path: Path | None, token_file: Path | None) -> None:
    """Prompt explicitly and create a private token file without displaying it."""
    config = scrontab._read_config(config_path) if config_path is not None else {}
    path = token_file or Path(config.get("E3SM_DIAGS_TOKEN_FILE", ""))
    if not path.is_absolute() or str(path).startswith("/absolute/path"):
        raise ValueError(
            "Set E3SM_DIAGS_TOKEN_FILE in configuration or specify TOKEN_FILE=/path."
        )
    if not sys.stdin.isatty():
        raise ValueError(
            "Token creation is interactive; run make ops-token-create from a terminal."
        )
    if path.exists() or path.is_symlink():
        raise FileExistsError(
            f"Token file already exists: {path}. For rotation, create a new TOKEN_FILE "
            "and update E3SM_DIAGS_TOKEN_FILE in controller.env after review."
        )
    token = getpass.getpass("Paste the E3SM Diags token: ")
    if not token:
        raise ValueError("Token must not be empty.")
    path.parent.mkdir(mode=0o700, parents=True, exist_ok=True)
    descriptor = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_EXCL, 0o600)
    with os.fdopen(descriptor, "w", encoding="utf-8") as stream:
        stream.write(token + "\n")
    _output(f"Token created: {path} (mode 0600)")


def _dispatch(args: argparse.Namespace) -> None:
    """Dispatch consolidated operations actions; require confirmation first."""
    if args.command == "help":
        _output(_HELP)
        return
    if args.command in ("run", "report", "enable", "disable") or (
        args.command == "env" and args.action == "update"
    ):
        if args.confirm != "YES":
            raise ValueError(f"Refusing {args.command}; specify CONFIRM=YES.")
    if args.command == "init":
        if args.operations_dir is None or not args.operations_dir.is_absolute():
            raise ValueError("Initial setup requires OPERATIONS_DIR=/absolute/path.")
        checkout, created_config = scrontab.initialize_operations(
            args.operations_dir, args.repository_url, args.branch
        )
        _output(
            f"Checkout: {checkout}\nConfiguration: {created_config}\nEdit configuration, then run ops-env ACTION=create and ops-enable CONFIRM=YES. Scheduling is not enabled."
        )
        return
    if args.command == "disable":
        scrontab.remove_scrontab()
        return
    if args.command == "token-create" and args.token_file is not None:
        _token_create(None, args.token_file)
        return
    config_path = resolve_config(args.config)
    if args.command == "status":
        dashboard(config_path)
    elif args.command == "logs":
        logs(config_path, job=args.job, lines=args.lines)
    elif args.command == "update":
        update_checkout(config_path)
    elif args.command in ("run", "report"):
        config = scrontab._read_config(config_path)
        component = "controller" if args.command == "run" else "reporter"
        script = (
            Path(config["REPOSITORY"])
            / "tests"
            / "complete_run"
            / f"complete-run-{component}.sh"
        )
        subprocess.run([str(script), str(config_path), "--manual"], check=True)
    elif args.command == "env":
        if args.action == "create":
            scrontab.create_controller_environment(config_path)
        elif args.action == "update":
            _output(
                f"Provenance: {scrontab.update_controller_environment(config_path, confirmed=True)}"
            )
        else:
            raise ValueError(
                "Specify ACTION=create or ACTION=update; environment actions are never automatic."
            )
    elif args.command == "enable":
        environment = scrontab._controller_environment(config_path)
        with scrontab._controller_lock(environment.results_root):
            validate_deployment(config_path)
            scrontab.install_scrontab(config_path)
    elif args.command == "token-create":
        _token_create(config_path, args.token_file)
    elif args.command == "shortcut":
        _output(shortcut(config_path))


def main(argv: Sequence[str] | None = None) -> int:
    """Run the consolidated operator CLI with actionable diagnostics."""
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog=_HELP,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("command", choices=(*_EVERYDAY, *_ADMIN))
    parser.add_argument(
        "--config", type=Path, default=os.environ.get("OPS_INPUT_CONFIG") or None
    )
    parser.add_argument("--confirm", default=os.environ.get("OPS_INPUT_CONFIRM", ""))
    parser.add_argument("--action", default=os.environ.get("OPS_INPUT_ACTION", ""))
    parser.add_argument("--job", default=os.environ.get("OPS_INPUT_JOB") or None)
    parser.add_argument(
        "--lines", type=int, default=os.environ.get("OPS_INPUT_LINES") or "100"
    )
    parser.add_argument(
        "--operations-dir",
        type=Path,
        default=os.environ.get("OPS_INPUT_OPERATIONS_DIR") or None,
    )
    parser.add_argument(
        "--repository-url",
        default=os.environ.get("OPS_INPUT_REPOSITORY_URL")
        or "https://github.com/E3SM-Project/e3sm_diags.git",
    )
    parser.add_argument(
        "--branch", default=os.environ.get("OPS_INPUT_BRANCH") or "main"
    )
    parser.add_argument(
        "--token-file",
        type=Path,
        default=os.environ.get("OPS_INPUT_TOKEN_FILE") or None,
    )
    args = parser.parse_args(argv)
    try:
        _dispatch(args)
    except (
        OSError,
        ValueError,
        RuntimeError,
        KeyError,
        subprocess.CalledProcessError,
    ) as error:
        detail = (
            error.stderr if isinstance(error, subprocess.CalledProcessError) else None
        )
        parser.error((detail or str(error)).strip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
