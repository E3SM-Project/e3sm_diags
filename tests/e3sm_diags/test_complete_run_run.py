"""Tests for complete-run result provenance."""

from __future__ import annotations

import subprocess
from argparse import Namespace
from pathlib import Path

import pytest

from tests.complete_run import run


def test_export_environment_writes_conda_output(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    monkeypatch.setattr(
        run.subprocess,
        "run",
        lambda *_, **__: subprocess.CompletedProcess([], 0, "name: test\n"),
    )

    path = run._export_environment(results_dir)

    assert path == results_dir / "prov" / "environment.yml"
    assert path.read_text(encoding="utf-8") == "name: test\n"


def test_export_environment_failure_does_not_publish_manifest_provenance(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    results_dir = tmp_path / "results"
    results_dir.mkdir()
    monkeypatch.setattr(
        run.subprocess,
        "run",
        lambda *_, **__: (_ for _ in ()).throw(
            subprocess.CalledProcessError(1, "conda")
        ),
    )

    with pytest.raises(RuntimeError, match="Unable to export"):
        run._export_environment(results_dir)

    assert not (results_dir / "prov" / "environment.yml").exists()


@pytest.fixture
def complete_run_args(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Namespace:
    args = run._build_parser().parse_args(["--results-dir", str(tmp_path / "results")])
    monkeypatch.setattr(run, "_validate_input_paths", lambda _: None)
    monkeypatch.setattr(run, "build_complete_run_params", lambda _: [])
    monkeypatch.setattr(run, "_build_manifest", lambda *_, **__: {})
    monkeypatch.setattr(run, "make_tree_public", lambda _: None)
    return args


@pytest.mark.parametrize("driver_exports", [True, False])
def test_run_finalizes_environment_provenance(
    complete_run_args: Namespace,
    monkeypatch: pytest.MonkeyPatch,
    driver_exports: bool,
) -> None:
    results_dir = Path(complete_run_args.results_dir)
    environment_path = results_dir / "prov" / "environment.yml"
    exports: list[Path] = []

    def fake_runner(_: object) -> list:
        environment_path.parent.mkdir(parents=True)
        if driver_exports:
            environment_path.write_text("name: driver\n", encoding="utf-8")
        return []

    def fake_export(path: Path) -> Path:
        exports.append(path)
        environment_path.write_text("name: fallback\n", encoding="utf-8")
        return environment_path

    monkeypatch.setattr(run.runner, "run_diags", fake_runner)
    monkeypatch.setattr(run, "_export_environment", fake_export)

    assert run._run_complete_run(complete_run_args) == []
    assert exports == ([] if driver_exports else [results_dir])
    assert environment_path.read_text(encoding="utf-8") == (
        "name: driver\n" if driver_exports else "name: fallback\n"
    )
    assert (results_dir / run._MANIFEST_FILENAME).is_file()


@pytest.mark.parametrize(
    ("preexisting", "kind"),
    [
        (True, "file"),
        (True, "empty"),
        (True, "directory"),
        (True, "symlink"),
        (False, "empty"),
        (False, "directory"),
        (False, "symlink"),
    ],
)
def test_run_rejects_preexisting_or_invalid_environment(
    complete_run_args: Namespace,
    monkeypatch: pytest.MonkeyPatch,
    kind: str,
    preexisting: bool,
) -> None:
    results_dir = Path(complete_run_args.results_dir)
    environment_path = results_dir / "prov" / "environment.yml"

    def create_provenance() -> None:
        environment_path.parent.mkdir(parents=True)
        if kind == "directory":
            environment_path.mkdir()
        elif kind == "symlink":
            environment_path.symlink_to(results_dir / "missing.yml")
        else:
            environment_path.write_text(
                "name: existing\n" if kind == "file" else "", encoding="utf-8"
            )

    def fake_runner(_: object) -> list:
        if preexisting:
            pytest.fail("Diagnostics must not run with preexisting provenance")
        create_provenance()
        return []

    if preexisting:
        create_provenance()
    monkeypatch.setattr(run.runner, "run_diags", fake_runner)
    monkeypatch.setattr(
        run,
        "_export_environment",
        lambda _: pytest.fail("Existing provenance must not be overwritten"),
    )

    expected_error = FileExistsError if preexisting else RuntimeError
    with pytest.raises(expected_error, match="environment provenance"):
        run._run_complete_run(complete_run_args)

    assert not (results_dir / run._MANIFEST_FILENAME).exists()
    if kind == "file":
        assert environment_path.read_text(encoding="utf-8") == "name: existing\n"
