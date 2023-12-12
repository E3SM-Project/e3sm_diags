import os
import re
import shutil
import sys
from typing import List

import pytest
from PIL import Image, ImageChops, ImageDraw

from e3sm_diags.logger import custom_logger
from e3sm_diags.run import runner
from tests.integration.config import TEST_IMAGES_PATH, TEST_ROOT_PATH
from tests.integration.utils import _get_test_params

CFG_PATH = os.path.join(TEST_ROOT_PATH, "all_sets.cfg")


logger = custom_logger(__name__)


@pytest.fixture(scope="module")
def run_diags_and_get_results_dir() -> str:
    """Run the diagnostics and get the results directory containing the images.

    The scope of this fixture is at the module level so that it only runs
    once, then each individual test can reference the result directory.

    Returns
    -------
    str
        The path to the results directory.
    """
    # Set -d flag to use the .cfg file for running additional diagnostic sets.
    sys.argv.extend(["-d", CFG_PATH])

    params = _get_test_params()
    results = runner.run_diags(params)

    results_dir = results[0].results_dir

    logger.info(f"results_dir={results_dir}")

    return results_dir


class TestAllSetsImageDiffs:
    @pytest.fixture(autouse=True)
    def setup(self, run_diags_and_get_results_dir):
        self.results_dir = run_diags_and_get_results_dir

    def test_num_images_is_the_same_as_the_expected(self):
        actual_num_images, actual_images = self._count_images_in_dir(
            f"{TEST_ROOT_PATH}/all_sets_results_test"
        )
        expected_num_images, expected_images = self._count_images_in_dir(
            TEST_IMAGES_PATH
        )

        assert actual_images == expected_images
        assert actual_num_images == expected_num_images

    def test_area_mean_time_series_plot_diffs(self):
        set_name = "area_mean_time_series"
        variables = ["TREFHT"]
        for variable in variables:
            variable_lower = variable.lower()

            # Check PNG path is the same as the expected.
            png_path = "{}/{}.png".format(set_name, variable)
            full_png_path = os.path.join(self.results_dir, png_path)
            path_exists = os.path.exists(full_png_path)

            assert path_exists

            # Check full HTML path is the same as the expected.
            filename = "viewer/{}/variable/{}/plot.html".format(
                set_name, variable_lower
            )
            html_path = os.path.join(self.results_dir, filename)
            self._check_html_image(html_path, png_path, full_png_path)

    def test_cosp_histogram_plot_diffs(self):
        self._check_plots_generic(
            set_name="cosp_histogram",
            case_id="MISR-COSP",
            ref_name="MISRCOSP",
            variables=["COSP_HISTOGRAM_MISR"],
            region="global",
        )

    def test_enso_diags_map_diffs(self):
        case_id = "TREFHT-response-map"
        self._check_enso_map_plots(case_id)

    def test_enso_diags_map_with_start_yrs_diffs(self):
        case_id = "TREFHT-response-map-start-yrs"
        self._check_enso_map_plots(case_id)

    def test_enso_diags_map_test_with_ref_yrs_diffs(self):
        case_id = "TREFHT-response-map-test-ref-yrs"
        self._check_enso_map_plots(case_id)

    def test_enso_diags_scatter_plot_diffs(self):
        case_id = "TREFHT-response-scatter"
        self._check_enso_scatter_plots(case_id)

    def test_enso_diags_scatter_with_start_yrs_plot_diffs(self):
        case_id = "TREFHT-response-scatter-start-yrs"
        self._check_enso_scatter_plots(case_id)

    def test_enso_diags_scatter_with_test_ref_yrs_plot_diffs(self):
        case_id = "TREFHT-response-scatter-test-ref-yrs"
        self._check_enso_scatter_plots(case_id)

    def test_lat_lon_plot_diffs(self):
        self._check_plots_plevs("lat_lon", "global", [850.0])

    def test_lat_lon_regional_plot_diffs(self):
        self._check_plots_plevs("lat_lon", "CONUS_RRM", [850.0])

    def test_meridional_mean_2d_plot_diffs(self):
        self._check_plots_2d("meridional_mean_2d")

    def test_polar_plot_diffs(self):
        self._check_plots_plevs("polar", "polar_S", [850.0])

    def test_qbo_plot_diffs(self):
        case_id = "qbo-test"
        case_id_lower = case_id.lower()
        set_name = "qbo"

        # Check PNG path is the same as the expected.
        png_path = "{}/{}/qbo_diags.png".format(set_name, case_id)
        full_png_path = os.path.join(self.results_dir, png_path)
        path_exists = os.path.exists(full_png_path)

        assert path_exists

        # Check full HTML path is the same as the expected.
        # viewer/qbo/variable/era-interim/plot.html
        filename = "viewer/{}/variable/{}/plot.html".format(set_name, case_id_lower)
        html_path = os.path.join(self.results_dir, filename)
        self._check_html_image(html_path, png_path, full_png_path)

    def test_streamflow_plot_diffs(self):
        self._check_streamflow_plots()

    def test_zonal_mean_2d_plot_diffs(self):
        self._check_plots_2d("zonal_mean_2d")

    def test_zonal_mean_xy_plot_diffs(self):
        self._check_plots_plevs("zonal_mean_xy", "global", [200.0])

    # Utility methods.
    # --------------------------------------------------------------------------
    def _count_images_in_dir(self, directory):
        images = []
        for root, _, filenames in os.walk(directory):
            # download_data.py won't download files in the viewer directory
            # because the webpage is more than a simple page of links.
            if "viewer" not in root:
                for file in filenames:
                    if file.endswith(".png"):
                        images.append(file)
        return len(images), images

    def _check_html_image(self, html_path, png_path, full_png_path):
        # Check HTML image tags exist.
        img_src = None
        option_value = None
        href = None
        with open(html_path, "r") as html:
            for line in html:
                # If `img_src` is not defined yet:
                if not img_src:
                    re_str = '<img src="../../../../{}">'.format(png_path)
                    img_src = re.search(re_str, line)
                # If `option_value` is not defined yet:
                if not option_value:
                    re_str = '<option value="../../../../{}">'.format(png_path)
                    option_value = re.search(re_str, line)
                # If `href` is not defined yet:
                if not href:
                    re_str = 'href="../../../../{}">'.format(png_path)
                    href = re.search(re_str, line)

        assert img_src is not None
        assert option_value is not None
        assert href is not None

        image_name = os.path.split(png_path)[-1]
        path_to_actual_png = full_png_path
        path_to_expected_png = f"{TEST_IMAGES_PATH}/{png_path}"

        if "CHECK_IMAGES" in os.environ:
            # Set `export CHECK_IMAGES=True` to do a pixel-by-pixel image comparison check.
            check_images = os.environ["CHECK_IMAGES"].lower() in [
                "true",
                "yes",
                "t",
                "y",
                "1",
            ]
        else:
            check_images = True

        if check_images:
            mismatched_images: List[str] = []
            _compare_images(
                mismatched_images,
                image_name,
                path_to_actual_png,
                path_to_expected_png,
            )
            assert len(mismatched_images) == 0

    def _check_plots_generic(
        self, set_name, case_id, ref_name, variables, region, plev=None
    ):
        case_id_lower = case_id.lower()
        ref_name_lower = ref_name.lower()
        region_lower = region.lower()
        seasons = ["ANN"]
        for variable in variables:
            variable_lower = variable.lower()
            for season in seasons:
                season_lower = season.lower()

                # Check PNG path is the same as the expected.

                if plev:
                    # 200.9 would just show up as 200 in the file paths.
                    plev_str = "%.0f" % plev
                    plev_png_str = "{}-".format(plev_str)
                    plev_html_str = "{}mb-".format(plev_str)
                else:
                    plev_png_str = ""
                    plev_html_str = ""
                png_path = "{}/{}/{}-{}-{}{}-{}.png".format(
                    set_name,
                    case_id,
                    ref_name,
                    variable,
                    plev_png_str,
                    season,
                    region,
                )

                full_png_path = os.path.join(self.results_dir, png_path)
                path_exists = os.path.exists(full_png_path)

                assert path_exists

                # Check full HTML path is the same as the expected.
                filename = "viewer/{}/{}/{}-{}{}-{}/{}.html".format(
                    set_name,
                    case_id_lower,
                    variable_lower,
                    plev_html_str,
                    region_lower,
                    ref_name_lower,
                    season_lower,
                )
                html_path = os.path.join(self.results_dir, filename)
                self._check_html_image(html_path, png_path, full_png_path)

    def _check_plots_2d(self, set_name):
        self._check_plots_generic(
            set_name=set_name,
            case_id="ERA-Interim",
            ref_name="ERA-Interim",
            variables=["T"],
            region="global",
        )

    def _check_plots_plevs(self, set_name, region, plevs):
        for plev in plevs:
            self._check_plots_generic(
                set_name=set_name,
                case_id="ERA-Interim",
                ref_name="ERA-Interim",
                variables=["T"],
                region=region,
                plev=plev,
            )

    def _check_enso_map_plots(self, case_id):
        case_id_lower = case_id.lower()
        nino_region_lower = "NINO34".lower()
        set_name = "enso_diags"
        variables = ["TREFHT"]

        for variable in variables:
            variable_lower = variable.lower()

            # Check PNG path is the same as the expected.
            png_path = "{}/{}/regression-coefficient-{}-over-{}.png".format(
                set_name, case_id, variable_lower, nino_region_lower
            )
            full_png_path = os.path.join(self.results_dir, png_path)
            path_exists = os.path.exists(full_png_path)

            assert path_exists

            # Check full HTML path is the same as the expected.
            filename = "viewer/{}/map/{}/plot.html".format(set_name, case_id_lower)
            html_path = os.path.join(self.results_dir, filename)
            self._check_html_image(html_path, png_path, full_png_path)

    def _check_enso_scatter_plots(self, case_id):
        case_id_lower = case_id.lower()
        set_name = "enso_diags"
        variables = ["TREFHT"]

        for variable in variables:
            region = "NINO3"

            # Check PNG path is the same as the expected.
            png_path = "{}/{}/feedback-{}-{}-TS-NINO3.png".format(
                set_name, case_id, variable, region
            )
            full_png_path = os.path.join(self.results_dir, png_path)
            path_exists = os.path.exists(full_png_path)

            assert path_exists

            # Check full HTML path is the same as the expected.
            filename = "viewer/{}/scatter/{}/plot.html".format(set_name, case_id_lower)
            html_path = os.path.join(self.results_dir, filename)
            self._check_html_image(html_path, png_path, full_png_path)

    def _check_streamflow_plots(self):
        case_id = "RIVER_DISCHARGE_OVER_LAND_LIQ_GSIM"
        case_id_lower = case_id.lower()
        set_name = "streamflow"
        variables = ["RIVER_DISCHARGE_OVER_LAND_LIQ"]
        for variable in variables:
            for plot_type in [
                "seasonality_map",
                "annual_map",
                "annual_scatter",
            ]:
                # Check PNG path is the same as the expected.
                png_path = "{}/{}/{}.png".format(set_name, case_id, plot_type)
                expected = (
                    "streamflow/RIVER_DISCHARGE_OVER_LAND_LIQ_GSIM/{}.png".format(
                        plot_type
                    )
                )
                assert png_path == expected

                # Check path exists
                full_png_path = os.path.join(self.results_dir, png_path)
                path_exists = os.path.exists(full_png_path)

                assert path_exists

                # Check HTML path is the same as the expected.
                if plot_type == "seasonality_map":
                    plot_label = "seasonality-map"
                elif plot_type == "annual_map":
                    plot_label = "mean-annual-streamflow-map"
                elif plot_type == "annual_scatter":
                    plot_label = "mean-annual-streamflow-scatter-plot"
                else:
                    raise Exception("Invalid plot_type={}".format(plot_type))
                html_path = "viewer/{}/{}/{}-{}/plot.html".format(
                    set_name, plot_label, case_id_lower, plot_type
                )
                expected = "viewer/streamflow/{}/river_discharge_over_land_liq_gsim-{}/plot.html".format(
                    plot_label, plot_type
                )
                assert html_path == expected

                # Check the full HTML path is the same as the expected.
                full_html_path = os.path.join(self.results_dir, html_path)
                self._check_html_image(full_html_path, png_path, full_png_path)


def _compare_images(
    mismatched_images: List[str],
    image_name: str,
    path_to_actual_png: str,
    path_to_expected_png: str,
) -> List[str]:
    # https://stackoverflow.com/questions/35176639/compare-images-python-pil

    actual_png = Image.open(path_to_actual_png).convert("RGB")
    expected_png = Image.open(path_to_expected_png).convert("RGB")
    diff = ImageChops.difference(actual_png, expected_png)

    diff_dir = f"{TEST_ROOT_PATH}image_check_failures"
    if not os.path.isdir(diff_dir):
        os.mkdir(diff_dir)

    bbox = diff.getbbox()
    # If `diff.getbbox()` is None, then the images are in theory equal
    if bbox is None:
        pass
    else:
        # Sometimes, a few pixels will differ, but the two images appear identical.
        # https://codereview.stackexchange.com/questions/55902/fastest-way-to-count-non-zero-pixels-using-python-and-pillow
        nonzero_pixels = (
            diff.crop(bbox)
            .point(lambda x: 255 if x else 0)
            .convert("L")
            .point(bool)
            .getdata()
        )
        num_nonzero_pixels = sum(nonzero_pixels)
        logger.info("\npath_to_actual_png={}".format(path_to_actual_png))
        logger.info("path_to_expected_png={}".format(path_to_expected_png))
        logger.info("diff has {} nonzero pixels.".format(num_nonzero_pixels))
        width, height = expected_png.size
        num_pixels = width * height
        logger.info("total number of pixels={}".format(num_pixels))
        fraction = num_nonzero_pixels / num_pixels
        logger.info("num_nonzero_pixels/num_pixels fraction={}".format(fraction))

        # Fraction of mismatched pixels should be less than 0.02%
        if fraction >= 0.0002:
            mismatched_images.append(image_name)

            simple_image_name = image_name.split("/")[-1].split(".")[0]
            shutil.copy(
                path_to_actual_png,
                os.path.join(diff_dir, "{}_actual.png".format(simple_image_name)),
            )
            shutil.copy(
                path_to_expected_png,
                os.path.join(diff_dir, "{}_expected.png".format(simple_image_name)),
            )
            # https://stackoverflow.com/questions/41405632/draw-a-rectangle-and-a-text-in-it-using-pil
            draw = ImageDraw.Draw(diff)
            (left, upper, right, lower) = diff.getbbox()
            draw.rectangle(((left, upper), (right, lower)), outline="red")
            diff.save(
                os.path.join(diff_dir, "{}_diff.png".format(simple_image_name)),
                "PNG",
            )

    return mismatched_images
