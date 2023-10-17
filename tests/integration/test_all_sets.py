import os
import re
import shutil
from typing import List

from tests.integration.config import TEST_DATA_DIR
from tests.integration.utils import run_cmd_and_pipe_stderr

# The path to the integration test data, which needs to be downloaded
# prior to running this test file.
MODULE_PATH = os.path.dirname(__file__)
TEST_DATA_PATH = os.path.join(MODULE_PATH, TEST_DATA_DIR)

# The path to the integration test diagnostics .cfg file.
CFG_PATH = os.path.join(MODULE_PATH, "all_sets_modified.cfg")
CFG_PATH = os.path.abspath(CFG_PATH)


class TestAllSets:
    def test_all_sets(self):
        expected_num_diags = 12

        # *_data_path needs to be added b/c the tests runs the diags from a different location
        # TODO: This test should be calling the Python modules directly rather
        # than passing the command to `run_cmd_and_pipe_stderr`. This allows
        # us to step directly in the module that is being called for debugging.
        # CoreParameter objects should be passed to the Run object.
        cmd = (
            f"e3sm_diags_driver.py -d {CFG_PATH} "
            f"--reference_data_path {TEST_DATA_PATH} "
            f"--test_data_path {TEST_DATA_PATH}"
        )

        stderr = run_cmd_and_pipe_stderr(cmd)

        # count the number of pngs in viewer_dir
        results_dir = self._get_results_dir(stderr)
        count = self._count_images(results_dir)

        # -1 is needed because of the E3SM logo in the viewer html
        assert count - 1 == expected_num_diags

        shutil.rmtree(results_dir)  # remove all generated results from the diags

    def _get_results_dir(self, output: List[str]):
        """Given output from e3sm_diags_driver, extract the path to results_dir."""
        for line in output:
            match = re.search("Viewer HTML generated at (.*)viewer.*.html", line)
            if match:
                results_dir = match.group(1)
                return results_dir

        raise RuntimeError("No viewer directory listed in output: {}".format(output))

    def _count_images(self, directory: str):
        """Count the number of images of type file_type in directory"""
        count = 0
        for _, __, files in os.walk(directory):
            for f in files:
                if f.endswith("png"):
                    count += 1
        return count
