"""The complete E3SM Diagnostic run.

Due to the large amount of data required to run, this test will be run manually
on Anvil (rather than as part of the CI tests).

Run the following first:
  - srun --pty --nodes=1 --time=01:00:00 /bin/bash
  - source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
  - Or: source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh
"""
import os

# This test should be run with the latest E3SM Diags tutorial code.
from examples.run_v2_6_0_all_sets_E3SM_machines import run_lcrc
from tests.integration.utils import _compare_images


class TestCompleteRun:
    def test_complete_run(self):
        actual_images_dir = run_lcrc(".")

        # The expected_images_file lists all images we expect to compare.
        expected_images_file = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/image_list_all_sets.txt"
        expected_images_dir = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/all_sets"

        mismatched_images = []  # type:ignore

        with open(expected_images_file) as f:
            for line in f:
                image_name = line.strip("./").strip("\n")
                path_to_actual_png = os.path.join(actual_images_dir, image_name)
                path_to_expected_png = os.path.join(expected_images_dir, image_name)

                mismatched_images = _compare_images(
                    mismatched_images,
                    image_name,
                    path_to_actual_png,
                    path_to_expected_png,
                )

        assert len(mismatched_images) == 0
