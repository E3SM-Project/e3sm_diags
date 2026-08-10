from __future__ import annotations

from pathlib import Path

from PIL import Image

from tests.integration import refresh_plot_image_baselines as refresh_baselines
from tests.integration.plot_image_regression_case import ImageRegressionCase


def _write_png(path: Path, differing_pixels: int = 0) -> None:
    image = Image.new("RGB", (10, 10), "black")

    for pixel_index in range(differing_pixels):
        image.putpixel((pixel_index, 0), (255, 255, 255))

    image.save(path)


def _case(
    baseline_dir: Path, source_path: Path, threshold: float
) -> ImageRegressionCase:
    return ImageRegressionCase(
        case_id="test",
        baseline_dir=baseline_dir,
        baseline_metadata_path=baseline_dir / "baseline_metadata.json",
        expected_image_filenames=("plot.png",),
        render=lambda _output_dir: (source_path,),
        mismatch_threshold=threshold,
    )


def test_refresh_keeps_images_within_mismatch_threshold(
    tmp_path: Path, monkeypatch
) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    expected_path = baseline_dir / "plot.png"
    source_path = tmp_path / "actual.png"
    _write_png(expected_path)
    _write_png(source_path, differing_pixels=1)
    metadata_writes: list[Path] = []
    monkeypatch.setattr(
        refresh_baselines,
        "write_runtime_metadata",
        lambda path: metadata_writes.append(Path(path)),
    )

    refresh_baselines.refresh_case_baselines(
        _case(baseline_dir, source_path, threshold=0.02)
    )

    assert Image.open(expected_path).getpixel((0, 0)) == (0, 0, 0)
    assert metadata_writes == []


def test_refresh_replaces_images_over_mismatch_threshold(
    tmp_path: Path, monkeypatch
) -> None:
    baseline_dir = tmp_path / "baselines"
    baseline_dir.mkdir()
    expected_path = baseline_dir / "plot.png"
    source_path = tmp_path / "actual.png"
    _write_png(expected_path)
    _write_png(source_path, differing_pixels=3)
    metadata_writes: list[Path] = []
    monkeypatch.setattr(
        refresh_baselines,
        "write_runtime_metadata",
        lambda path: metadata_writes.append(Path(path)),
    )

    refresh_baselines.refresh_case_baselines(
        _case(baseline_dir, source_path, threshold=0.02)
    )

    assert Image.open(expected_path).getpixel((0, 0)) == (255, 255, 255)
    assert metadata_writes == [baseline_dir / "baseline_metadata.json"]
