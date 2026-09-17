"""Tests for complete-run result provenance."""

from __future__ import annotations

import subprocess
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
