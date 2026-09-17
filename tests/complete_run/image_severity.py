"""Classify complete-run PNG differences by review severity.

Raw pixel counts are sensitive to anti-aliasing and small text-metric changes
from rendering-library upgrades. This module tolerates small local movements
while retaining changes to plotted content, layout, and compact text such as a
printed statistic.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Final, Literal

import numpy as np
from PIL import Image
from scipy import ndimage

ImageSeverity = Literal[
    "IDENTICAL", "NEGLIGIBLE", "MINOR", "MODERATE", "MAJOR", "STRUCTURAL"
]

IDENTICAL: Final[ImageSeverity] = "IDENTICAL"
NEGLIGIBLE: Final[ImageSeverity] = "NEGLIGIBLE"
MINOR: Final[ImageSeverity] = "MINOR"
MODERATE: Final[ImageSeverity] = "MODERATE"
MAJOR: Final[ImageSeverity] = "MAJOR"
STRUCTURAL: Final[ImageSeverity] = "STRUCTURAL"

SEVERITY_ORDER: Final[tuple[ImageSeverity, ...]] = (
    IDENTICAL,
    NEGLIGIBLE,
    MINOR,
    MODERATE,
    MAJOR,
    STRUCTURAL,
)
REVIEWABLE_SEVERITIES: Final[tuple[ImageSeverity, ...]] = (
    MINOR,
    MODERATE,
    MAJOR,
    STRUCTURAL,
)

# These values are the zppy #865 calibration, which was validated against
# thousands of e3sm_diags images. They must be re-evaluated on complete-run
# baselines before being relied upon for an environment-regression gate.
# Source: https://github.com/E3SM-Project/zppy/pull/865
SHIFT_TOLERANCE_PIXELS: Final = 4
INTENSITY_TOLERANCE: Final = 32
NEGLIGIBLE_MAX: Final = 0.005
MINOR_MAX: Final = 0.02
MODERATE_MAX: Final = 0.08
STRUCTURAL_GEOMETRY_CHANGE: Final = 0.10
NOTABLE_GEOMETRY_CHANGE: Final = 0.006
BACKGROUND_TOLERANCE: Final = 6
LOCALIZED_SHIFT_TOLERANCE_PIXELS: Final = 1
LOCALIZED_INTENSITY_TOLERANCE: Final = 80
LOCALIZED_MIN_SPOT_PIXELS: Final = 8
LOCALIZED_MIN_SPOT_EXTENT: Final = 3
LOCALIZED_MIN_TOTAL_PIXELS: Final = 20
LOCALIZED_MAX_SPOTS: Final = 20


@dataclass(frozen=True)
class ImageComparison:
    """Severity result for one shared PNG pair."""

    relative_path: Path
    severity: ImageSeverity
    content_fraction: float
    geometry_change: float
    actual_size: tuple[int, int]
    baseline_size: tuple[int, int]
    cause: str
    localized_pixels: int = 0
    raw_fraction: float = 0.0

    @property
    def needs_review(self) -> bool:
        """Whether this result should fail comparison and get artifacts."""
        return self.severity in REVIEWABLE_SEVERITIES

    @property
    def is_cosmetic(self) -> bool:
        """Whether this result differs only by rendering noise."""
        return self.severity == NEGLIGIBLE

    def sort_key(self) -> tuple[int, float]:
        """Return a worst-first ordering key."""
        return (-SEVERITY_ORDER.index(self.severity), -self.content_fraction)


def compare_pngs(
    dev_path: str | Path,
    baseline_path: str | Path,
    *,
    relative_path: str | Path,
) -> ImageComparison:
    """Compare one PNG pair and classify its severity."""
    dev_file = Path(dev_path)
    baseline_file = Path(baseline_path)
    relative = Path(relative_path)

    if _files_are_identical(dev_file, baseline_file):
        size = _image_size(baseline_file)
        return ImageComparison(relative, IDENTICAL, 0.0, 0.0, size, size, "no change")

    actual = _load_rgb(dev_file)
    expected = _load_rgb(baseline_file)
    actual_size = actual.shape[:2]
    expected_size = expected.shape[:2]
    if actual.shape == expected.shape and np.array_equal(actual, expected):
        return ImageComparison(
            relative, IDENTICAL, 0.0, 0.0, actual_size, expected_size, "no change"
        )

    raw_fraction = _raw_difference_fraction(actual, expected)
    actual = trim_background(actual)
    expected = trim_background(expected)
    geometry_change = _relative_size_change(actual.shape[:2], expected.shape[:2])
    cause = _describe_cause(geometry_change, actual.shape[:2], expected.shape[:2])
    if geometry_change >= STRUCTURAL_GEOMETRY_CHANGE:
        return ImageComparison(
            relative,
            STRUCTURAL,
            1.0,
            geometry_change,
            actual_size,
            expected_size,
            cause,
            raw_fraction=raw_fraction,
        )

    localized_pixels = 0
    if actual.shape == expected.shape:
        localized_pixels = localized_change_pixels(actual, expected)
    if actual.shape != expected.shape:
        height, width = expected.shape[:2]
        actual = np.asarray(
            Image.fromarray(actual).resize((width, height), Image.Resampling.BILINEAR)
        )

    difference = tolerant_difference(actual, expected)
    content_fraction = float((difference > INTENSITY_TOLERANCE).mean())
    severity = _band(content_fraction)

    if geometry_change >= NOTABLE_GEOMETRY_CHANGE:
        severity = _promote(severity)
    if localized_pixels >= LOCALIZED_MIN_TOTAL_PIXELS and severity == NEGLIGIBLE:
        severity = MINOR
        cause = "small isolated difference"

    return ImageComparison(
        relative,
        severity,
        content_fraction,
        geometry_change,
        actual_size,
        expected_size,
        cause,
        localized_pixels,
        raw_fraction,
    )


def trim_background(image: np.ndarray) -> np.ndarray:
    """Trim a near-uniform border introduced by tight figure bounding boxes."""
    background = image[0, 0].astype(np.int16)
    is_content = (
        np.abs(image.astype(np.int16) - background).max(axis=2) > BACKGROUND_TOLERANCE
    )

    if not is_content.any():
        return image

    rows = np.where(is_content.any(axis=1))[0]
    columns = np.where(is_content.any(axis=0))[0]

    return image[rows[0] : rows[-1] + 1, columns[0] : columns[-1] + 1]


def tolerant_difference(
    actual: np.ndarray,
    expected: np.ndarray,
    radius: int = SHIFT_TOLERANCE_PIXELS,
) -> np.ndarray:
    """Return bidirectional per-pixel differences while forgiving local shifts."""
    size = 2 * radius + 1
    worst = np.zeros(actual.shape[:2], dtype=np.float32)
    for channel in range(3):
        actual_channel = actual[..., channel].astype(np.float32)
        expected_channel = expected[..., channel].astype(np.float32)
        actual_vs_expected = np.maximum(
            actual_channel - ndimage.grey_dilation(expected_channel, size=size),
            ndimage.grey_erosion(expected_channel, size=size) - actual_channel,
        )
        expected_vs_actual = np.maximum(
            expected_channel - ndimage.grey_dilation(actual_channel, size=size),
            ndimage.grey_erosion(actual_channel, size=size) - expected_channel,
        )
        worst = np.maximum(
            worst,
            np.maximum(actual_vs_expected, expected_vs_actual).clip(min=0),
        )
    return worst


def localized_change_pixels(actual: np.ndarray, expected: np.ndarray) -> int:
    """Return compact, high-contrast differences that can represent changed text."""
    difference = tolerant_difference(actual, expected, LOCALIZED_SHIFT_TOLERANCE_PIXELS)
    strong = difference > LOCALIZED_INTENSITY_TOLERANCE

    if not strong.any():
        return 0

    labels, count = ndimage.label(strong)
    if count == 0:
        return 0

    sizes = ndimage.sum(strong, labels, range(1, count + 1))
    boxes = ndimage.find_objects(labels)
    big_enough = [
        index for index, size in enumerate(sizes) if size >= LOCALIZED_MIN_SPOT_PIXELS
    ]
    if len(big_enough) > LOCALIZED_MAX_SPOTS:
        return 0

    total = 0
    for index in big_enough:
        box = boxes[index]

        if box is None:
            continue

        rows, columns = box
        if (
            min(rows.stop - rows.start, columns.stop - columns.start)
            < LOCALIZED_MIN_SPOT_EXTENT
        ):
            continue
        total += int(sizes[index])

    return total


def _files_are_identical(path_a: Path, path_b: Path) -> bool:
    if path_a.stat().st_size != path_b.stat().st_size:
        return False

    chunk_size = 1 << 16

    with path_a.open("rb") as file_a, path_b.open("rb") as file_b:
        while True:
            block_a = file_a.read(chunk_size)

            if block_a != file_b.read(chunk_size):
                return False
            if not block_a:
                return True


def _image_size(path: Path) -> tuple[int, int]:
    with Image.open(path) as image:
        return (image.height, image.width)


def _load_rgb(path: Path) -> np.ndarray:
    with Image.open(path) as image:
        return np.asarray(image.convert("RGB"))


def _raw_difference_fraction(actual: np.ndarray, expected: np.ndarray) -> float:
    height = min(actual.shape[0], expected.shape[0])
    width = min(actual.shape[1], expected.shape[1])

    if height == 0 or width == 0:
        return 1.0

    actual_shared = actual[:height, :width].astype(np.int16)
    expected_shared = expected[:height, :width].astype(np.int16)

    return float((np.abs(actual_shared - expected_shared).max(axis=2) > 0).mean())


def _relative_size_change(actual: tuple[int, int], expected: tuple[int, int]) -> float:
    return max(
        abs(actual[0] - expected[0]) / max(actual[0], expected[0], 1),
        abs(actual[1] - expected[1]) / max(actual[1], expected[1], 1),
    )


def _describe_cause(
    geometry_change: float,
    actual_size: tuple[int, int],
    expected_size: tuple[int, int],
) -> str:
    height_change = abs(expected_size[0] - actual_size[0])
    width_change = abs(expected_size[1] - actual_size[1])

    if geometry_change >= STRUCTURAL_GEOMETRY_CHANGE:
        return "figure size changed a lot"
    if height_change >= 20:
        return "figure height changed"
    if geometry_change >= NOTABLE_GEOMETRY_CHANGE or width_change > 3:
        return "figure size changed slightly"

    return "same size, content differs"


def _band(content_fraction: float) -> ImageSeverity:
    if content_fraction <= NEGLIGIBLE_MAX:
        return NEGLIGIBLE
    if content_fraction <= MINOR_MAX:
        return MINOR
    if content_fraction <= MODERATE_MAX:
        return MODERATE

    return MAJOR


def _promote(severity: ImageSeverity) -> ImageSeverity:
    index = SEVERITY_ORDER.index(severity)

    return SEVERITY_ORDER[min(index + 1, SEVERITY_ORDER.index(MAJOR))]
