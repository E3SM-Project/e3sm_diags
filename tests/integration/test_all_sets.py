import os
import re
import shutil
import unittest

from tests.integration.config import TEST_DATA_DIR
from tests.integration.utils import run_cmd_and_pipe_stderr


def count_images(directory, file_type="png"):
    """Count the number of images of type file_type in directory"""
    count = 0
    for _, __, files in os.walk(directory):
        for f in files:
            if f.endswith(file_type):
                count += 1
    return count


class TestAllSets(unittest.TestCase):
    def test_all_sets(self):
        pth = os.path.dirname(__file__)
        test_pth = os.path.join(pth, TEST_DATA_DIR)
        cfg_pth = os.path.join(pth, "all_sets_modified.cfg")
        cfg_pth = os.path.abspath(cfg_pth)
        expected_num_diags = 12

        # *_data_path needs to be added b/c the tests runs the diags from a different location
        # TODO: This test should be calling the Python modules directly rather
        # than passing the command to `run_cmd_and_pipe_stderr`. This allows
        # us to step directly in the module that is being called for debugging.
        # CoreParameter objects should be passed to the Run object.
        cmd = (
            f"e3sm_diags_driver.py -d {cfg_pth} "
            f"--reference_data_path {test_pth} "
            f"--test_data_path {test_pth}"
        )

        stderr = run_cmd_and_pipe_stderr(cmd)
        # count the number of pngs in viewer_dir
        results_dir = self._get_results_dir(stderr)
        count = count_images(results_dir)
        # -1 is needed because of the E3SM logo in the viewer html
        self.assertEqual(count - 1, expected_num_diags)

        shutil.rmtree(results_dir)  # remove all generated results from the diags

    def _get_results_dir(self, output):
        """Given output from e3sm_diags_driver, extract the path to results_dir."""
        for line in output:
            match = re.search("Viewer HTML generated at (.*)viewer.*.html", line)
            if match:
                results_dir = match.group(1)
                return results_dir

        self.fail("No viewer directory listed in output: {}".format(output))


if __name__ == "__main__":
    unittest.main()
