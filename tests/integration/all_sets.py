# Running the software with the API:
#    python all_sets.py -d all_sets.py
from e3sm_diags.run import Run
from tests.integration.utils import _get_test_params

run_object = Run()

params = _get_test_params()

run_object.run_diags(params)
