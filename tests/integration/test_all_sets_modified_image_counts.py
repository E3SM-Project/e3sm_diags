import os
import shutil
from typing import List

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner
from tests.integration.config import TEST_ROOT_PATH
from tests.integration.utils import _convert_cfg_to_param_objs, _count_images

# The path to the integration test diagnostics .cfg file.
CFG_PATH = os.path.join(TEST_ROOT_PATH, "all_sets_modified.cfg")
CFG_PATH = os.path.abspath(CFG_PATH)

# +1 is needed because of the E3SM logo that is used by the viewer HTML.
EXPECTED_NUM_IMAGES = 12 + 1


def test_all_sets_modified_produces_the_expected_number_of_images():
    params = _convert_cfg_to_param_objs(CFG_PATH)
    params_results: List[CoreParameter] = []

    for param in params:
        result = runner.run_diags([param], use_cfg=False)
        params_results.extend(result)

    # The result directory should be the same for all diagnostic sets.
    result_dir = params_results[0].results_dir
    num_images = _count_images(result_dir)

    assert num_images == EXPECTED_NUM_IMAGES

    # TODO: Result dir should be set thet temporary pytest location that
    # automatically gets cleaned up after every test.
    shutil.rmtree(result_dir)
