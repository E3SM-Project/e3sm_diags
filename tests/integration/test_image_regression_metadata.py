from __future__ import annotations

import json
from pathlib import Path

import pytest

from tests.integration import image_regression


@pytest.mark.parametrize("version", ["2026.09.1", None])
def test_write_runtime_metadata_records_uxarray_base(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch, version: str | None
) -> None:
    def get_version(dist_name: str) -> str:
        if dist_name == "uxarray" and version is not None:
            return version
        raise image_regression.metadata.PackageNotFoundError(dist_name)

    monkeypatch.setattr(image_regression.metadata, "version", get_version)
    monkeypatch.setattr(image_regression, "_MODULE_NAMES_BY_KEY", {})
    monkeypatch.setattr(image_regression, "_get_git_sha", lambda: "test-sha")

    output_path = tmp_path / image_regression.BASELINE_METADATA_FILENAME
    image_regression.write_runtime_metadata(output_path)

    recorded_metadata = json.loads(output_path.read_text(encoding="utf-8"))
    assert "uxarray" in recorded_metadata
    assert recorded_metadata["uxarray"] == version
