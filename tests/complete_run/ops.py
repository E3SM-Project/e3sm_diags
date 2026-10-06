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
from tests.complete_run.machine_paths import detect_machine_paths

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
  make ops-init [OPERATIONS_DIR=/path] [BRANCH=main] [REPOSITORY_URL=<url>]
  make ops-env ACTION=create
  make ops-env ACTION=update CONFIRM=YES
  make ops-enable CONFIRM=YES
  make ops-disable CONFIRM=YES
  make ops-token-create [TOKEN_FILE=/path]
  make ops-shortcut                 Print a shell function; do not install it

Configuration: CONFIG=/path/controller.env, then E3SM_DIAGS_OPS_CONFIG,
then an existing controller.env in this checkout's parent, then machine defaults.
On NERSC Perlmutter, operations and results paths are detected automatically.
Custom deployments can override CONFIG and OPERATIONS_DIR. Setup never enables scheduling.
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
    if not selected and not path.is_file():
        defaults = detect_machine_paths()
        if defaults is not None:
            path = defaults.operations_dir / "controller.env"
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(
            f"Controller configuration not found: {path}. Set CONFIG=/path/controller.env "
            "or E3SM_DIAGS_OPS_CONFIG; for initial setup use "
            "make ops-init (machine defaults) or make ops-init OPERATIONS_DIR=/path "
            "for a custom deployment, then edit controller.env."
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


def _section(title: str) -> None:
    """Separate dashboard sections without terminal-specific formatting."""
    _output(f"\n{title}\n{'-' * len(title)}")


def _field(label: str, value: object) -> None:
    """Print a labeled value, indenting multiline diagnostics."""
    lines = str(value if value is not None else "unknown").splitlines() or ["unknown"]
    _output(f"  {label}: {lines[0]}")
    for line in lines[1:]:
        _output(f"    {line}")


def _table(headers: tuple[str, ...], rows: Sequence[tuple[str, ...]]) -> None:
    """Render aligned columns with explicit spacing and no truncation."""
    widths = [max(len(row[i]) for row in [headers, *rows]) for i in range(len(headers))]
    for row in [headers, *rows]:
        _output(
            "  "
            + "  ".join(
                value.ljust(width) for value, width in zip(row, widths, strict=True)
            ).rstrip()
        )


def _component_name(command: str) -> str:
    """Use short names for known wrappers, preserving other job names."""
    name = Path(command).name
    for component in ("controller", "reporter"):
        if name == f"complete-run-{component}.sh":
            return component
    return command


def _installed_schedule() -> None:
    """Summarize only the installed managed block, never unrelated schedules."""
    _section("Installed schedule")
    _output("  Future recurring occurrences, not previous outcomes.")
    try:
        lines = scrontab._read_scrontab().splitlines()
    except (OSError, subprocess.CalledProcessError) as error:
        _field("Schedule unavailable", error)
        return
    begin, end = scrontab._MANAGED_BEGIN, scrontab._MANAGED_END
    if begin not in lines and end not in lines:
        _output("  No managed complete-run schedule installed.")
        return
    if (
        lines.count(begin) != 1
        or lines.count(end) != 1
        or lines.index(begin) >= lines.index(end)
    ):
        _output("  Invalid managed schedule markers; inspect with scrontab -l.")
        return
    rows = []
    for line in lines[lines.index(begin) + 1 : lines.index(end)]:
        entry = line.strip()
        if not entry or entry.startswith("#"):
            continue
        fields = entry.split(maxsplit=5)
        if len(fields) != 6:
            _field("Unrecognized schedule entry", entry)
            continue
        try:
            command = shlex.split(fields[5])[0]
        except (ValueError, IndexError):
            _field("Unrecognized schedule entry", entry)
            continue
        component = _component_name(command)
        if component not in ("controller", "reporter"):
            component = fields[5]
        expression = " ".join(fields[:5])
        when = {
            "0 13 * * 0": "Sunday 13:00",
            "0 14 * * 0": "Sunday 14:00",
            "0 16 * * 1": "Monday 16:00",
            "0 17 * * 1": "Monday 17:00",
        }.get(expression, "custom")
        rows.append((component, when, expression))
    if rows:
        _table(("Component", "When (UTC)", "Cron expression"), rows)
        _output(
            "  Known wrappers apply Pacific-time guards; controller also checks even ISO weeks."
        )
    else:
        _output("  No recognized cron entries in the managed block.")
    _output("  Full schedule and resource directives: scrontab -l")


def _cron_jobs() -> None:
    """Display Slurm eligible times with unambiguous column separators."""
    _section("Slurm cron jobs / next eligible occurrences")
    result = _inspect_command(
        [
            "squeue",
            "--me",
            "-q",
            "cron",
            "--noheader",
            "-O",
            "JobID:0|,Name:0|,State:0|,EligibleTime:0",
        ]
    )
    if result == "(none)":
        _output("  No cron jobs in the queue.")
        return
    if result.startswith("Unavailable:"):
        _output(f"  {result}")
        return
    rows = []
    for line in result.splitlines():
        fields = tuple(value.strip() for value in line.split("|"))
        if len(fields) != 4:
            _field("Unrecognized Slurm output", line)
            continue
        job, name, state, eligible = fields
        rows.append((job, _component_name(name), state, eligible))
    if rows:
        _table(("Job ID", "Component / name", "State", "Eligible time (Slurm)"), rows)
        _output(
            "  N/A means Slurm has no eligible time to display; it is not a run outcome."
        )


def dashboard(config_path: Path) -> None:
    """Show deployment health, schedule, jobs, and latest recorded run read-only."""
    try:
        config = scrontab._read_config(config_path)
    except (OSError, ValueError):
        # Validation below reports the diagnostic; schedule/queue inspection
        # remains useful even when deployment paths cannot be parsed.
        config = {}
    _output("E3SM Diagnostics operations (read-only)")
    _section("Deployment")
    _field("Configuration", config_path)
    for key, label in (
        ("REPOSITORY", "Repository"),
        ("LOG_DIR", "Logs"),
        ("CONDA_BASE", "Conda base"),
        ("CONTROLLER_ENV_PREFIX", "Controller environment"),
        ("RESULTS_ROOT", "Results root"),
    ):
        value = config.get(key, "")
        health = "ok" if value and Path(value).exists() else "missing"
        _field(label, f"[{health}] {value or '(not configured)'}")
    try:
        scrontab.validate_config(config_path)
        _field("Configuration health", "valid")
    except (OSError, ValueError) as error:
        _field("Configuration health", f"INVALID: {error}")
    repository = config.get("REPOSITORY")
    if repository:
        for component in ("controller", "reporter"):
            script = (
                Path(repository)
                / "tests"
                / "complete_run"
                / f"complete-run-{component}.sh"
            )
            executable = script.is_file() and os.access(script, os.X_OK)
            _field(
                f"{component.capitalize()} script",
                "[ok] executable"
                if executable
                else f"[missing or not executable] {script}",
            )
    _installed_schedule()
    _cron_jobs()
    _section("Latest automated run")
    try:
        _latest_run(config)
    except OSError as error:
        _field("Latest automated run unavailable", error)


def _latest_run(config: dict[str, str]) -> None:
    """Inspect only immediate automated run directories, never results trees."""
    root_value = config.get("RESULTS_ROOT")
    if not root_value:
        _field("Run", "none (RESULTS_ROOT is missing)")
        return
    root = Path(root_value) / "automation"
    candidates = sorted(
        (path for path in root.glob("*-????????-??????") if path.is_dir()),
        key=lambda path: (path.name.rsplit("-", 2)[-2:], path.name),
    )
    if not candidates:
        _field("Run", f"none in {root}")
        return
    run = candidates[-1]
    _field("Run", run.name)
    _field("Metadata directory", run)
    status = _metadata(run / "status.json")
    if "unavailable" in status:
        _field("Recorded outcome", status["unavailable"])
    else:
        _field("Recorded outcome", status.get("stage", "unknown"))
        _run_details(status)
    job = str(status.get("job_id", ""))
    if job.isdigit():
        _field("Recorded job ID", job)
        _field(
            "Recorded run accounting (not the next recurring occurrence)",
            _inspect_command(
                ["sacct", "-j", job, "--format=JobID,State,ExitCode", "--noheader"]
            ),
        )
    _publication_summary(run)


def _run_details(status: dict[str, object]) -> None:
    """Show useful run fields instead of an unbounded status JSON dump."""
    for key, label in (
        ("git_sha", "Commit"),
        ("submitted_at_utc", "Submitted (UTC)"),
        ("result_dir", "Results"),
        ("environment_name", "Environment"),
        ("environment_prefix", "Environment prefix"),
        ("error", "Error"),
    ):
        if status.get(key) is not None:
            _field(label, status[key])
    sets = status.get("selected_sets")
    if isinstance(sets, list):
        _field("Diagnostics", f"{len(sets)} selected sets (see status.json for names)")


def _publication_summary(run: Path) -> None:
    """Keep report, receipt, and failure states distinct, including diagnostics."""
    _section("Report and publication")
    report = _metadata(run / "automation-report.json")
    _field("Report", report.get("unavailable", report.get("status", "unknown")))
    publication = report.get("publication")
    if isinstance(publication, dict):
        _field("Publication", publication.get("status", "unknown"))
        if publication.get("discussion_url"):
            _field("Discussion", publication["discussion_url"])
    elif "unavailable" not in report:
        _field("Publication", "unknown (missing or invalid report publication field)")
    for filename, label in (
        ("publication-receipt.json", "Publication receipt"),
        ("publication-failure.json", "Publication failure"),
    ):
        if not (run / filename).is_file():
            continue
        payload = _metadata(run / filename)
        _field(label, payload.get("unavailable", payload.get("status", "unknown")))
        if payload.get("discussion_url"):
            _field(f"{label} Discussion", payload["discussion_url"])
    _output(
        "  Policy: qualifying failures only; a successful run need not have a Discussion."
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


def _operations_directory(explicit: Path | None) -> Path:
    """Select an absolute initialization directory without creating it.

    Raises
    ------
    ValueError
        If no machine default exists or the supplied path is relative.
    """
    directory = explicit
    if directory is None:
        defaults = detect_machine_paths()
        directory = defaults.operations_dir if defaults is not None else None
    if directory is None or not directory.is_absolute():
        raise ValueError(
            "No usable machine default or absolute operations directory; "
            "specify OPERATIONS_DIR=/absolute/path."
        )
    return directory


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
        operations_dir = _operations_directory(args.operations_dir)
        checkout, created_config = scrontab.initialize_operations(
            operations_dir, args.repository_url, args.branch
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
