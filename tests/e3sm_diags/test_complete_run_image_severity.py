"""Tests for severity scoring of complete-run PNG comparisons."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PIL import Image

from tests.complete_run.image_severity import (
    IDENTICAL,
    MAJOR,
    MINOR,
    NEGLIGIBLE,
    STRUCTURAL,
    ImageComparison,
    compare_pngs,
    tolerant_difference,
    trim_background,
)


def _blank(height: int = 80, width: int = 80) -> np.ndarray:
    return np.full((height, width, 3), 255, dtype=np.uint8)


def _write(path: Path, image: np.ndarray) -> Path:
    Image.fromarray(image).save(path)
    return path


class TestTolerantDifference:
    def test_forgives_a_small_shift(self):
        expected = _blank()
        expected[40, 10:70] = 0
        actual = _blank()
        actual[41, 10:70] = 0

        assert not (tolerant_difference(actual, expected) > 32).any()

    def test_detects_a_deleted_thin_feature_in_reverse_direction(self):
        expected = _blank()
        expected[40, 10:70] = 0

        assert (tolerant_difference(_blank(), expected) > 32).any()


class TestImageSeverity:
    def test_identical_image_is_distinct_from_cosmetic(self, tmp_path: Path):
        image = _blank()
        image[20:60, 20:60] = 0
        path = _write(tmp_path / "same.png", image)

        result = compare_pngs(path, path, relative_path="lat_lon/same.png")

        assert result.severity == IDENTICAL
        assert not result.needs_review
        assert not result.is_cosmetic

    def test_shifted_content_is_cosmetic(self, tmp_path: Path):
        expected = _blank(100, 100)
        expected[30:70, 30:70] = 0
        actual = _blank(100, 100)
        actual[31:71, 30:70] = 0

        result = compare_pngs(
            _write(tmp_path / "actual.png", actual),
            _write(tmp_path / "expected.png", expected),
            relative_path="lat_lon/shifted.png",
        )

        assert result.severity == NEGLIGIBLE
        assert result.is_cosmetic
        assert not result.needs_review

    def test_recolored_content_is_reviewable(self, tmp_path: Path):
        expected = _blank()
        expected[30:50, 30:50] = (255, 0, 0)
        actual = _blank()
        actual[30:50, 30:50] = (0, 0, 255)

        result = compare_pngs(
            _write(tmp_path / "actual.png", actual),
            _write(tmp_path / "expected.png", expected),
            relative_path="lat_lon/recolored.png",
        )

        # Background trimming compares the plot content, rather than its
        # whitespace margin. The entire remaining image was recolored, so it
        # is a MAJOR difference under the zppy-calibrated score bands.
        assert result.severity == MAJOR
        assert result.needs_review

    def test_large_layout_change_is_structural(self, tmp_path: Path):
        expected = _blank(100, 100)
        expected[20:80, 20:80] = 0
        actual = _blank(100, 100)
        actual[20:80, 20:60] = 0

        result = compare_pngs(
            _write(tmp_path / "actual.png", actual),
            _write(tmp_path / "expected.png", expected),
            relative_path="lat_lon/layout.png",
        )

        assert result.severity == STRUCTURAL
        assert result.needs_review

    def test_compact_changed_value_is_minor(self, tmp_path: Path):
        expected = _blank(200, 200)
        expected[40:160, 40:160] = 200
        actual = _blank(200, 200)
        actual[40:160, 40:160] = 200
        actual[100:106, 100:110] = 0

        result = compare_pngs(
            _write(tmp_path / "actual.png", actual),
            _write(tmp_path / "expected.png", expected),
            relative_path="lat_lon/statistic.png",
        )

        assert result.severity == MINOR
        assert result.needs_review

    def test_worst_first_sorting_prefers_severity_then_content(self):
        minor = ImageComparison(
            Path("minor.png"), MINOR, 0.01, 0.0, (1, 1), (1, 1), "content"
        )
        lower_major = ImageComparison(
            Path("lower.png"), MAJOR, 0.1, 0.0, (1, 1), (1, 1), "content"
        )
        higher_major = ImageComparison(
            Path("higher.png"), MAJOR, 0.2, 0.0, (1, 1), (1, 1), "content"
        )

        assert sorted(
            [minor, lower_major, higher_major], key=ImageComparison.sort_key
        ) == [
            higher_major,
            lower_major,
            minor,
        ]


def test_trim_background_removes_tight_bbox_margin():
    image = _blank(100, 101)
    image[40:60, 31:71] = 0

    assert trim_background(image).shape == (20, 40, 3)
