"""Run the 1044 spatial KeyError reproduction as a local e3sm_diags script.
Source: e3sm_diags lat_lon --no_viewer --reference_data_path '/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/' --test_data_path '/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr'  --results_dir '/lcrc/group/e3sm/public_html/diagnostic_output/ac.zhang40/tests/missing_bounds_TREFHT/latest_main' --case_id 'CRU_IPCC' --sets 'lat_lon' --variables 'TREFHT' --seasons 'ANN' --regions 'land_60S90N' --regrid_method 'bilinear' --multiprocessing --num_workers '8' --main_title 'TREFHT ANN land_60S90N' --contour_levels '-35' '-30' '-25' '-20' '-15' '-10' '-5' '0' '5' '10' '15' '20' '25' '30' '35' '40' --test_name 'v3.LR.historical_0051' --short_test_name 'v3.LR.historical_0051' --ref_name 'CRU' --reference_name 'CRU Global Monthly Mean T Land' --diff_title 'Model - Observations' --diff_levels '-15' '-10' '-5' '-2' '-1' '-0.5' '-0.2' '0.2' '0.5' '1' '2' '5' '10' '15'
"""

from __future__ import annotations

from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

logger = _setup_child_logger(__name__)


REFERENCE_DATA_PATH = "/lcrc/group/e3sm/diagnostics/observations/Atm/climatology/"
TEST_DATA_PATH = (
    "/lcrc/group/e3sm/ac.zhang40/zppy_example_v3.2.0/"
    "v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr"
)
RESULTS_DIR = (
    "/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/"
    "missing_bounds_TREFHT/latest_main"
)

CONTOUR_LEVELS = [
    -35.0,
    -30.0,
    -25.0,
    -20.0,
    -15.0,
    -10.0,
    -5.0,
    0.0,
    5.0,
    10.0,
    15.0,
    20.0,
    25.0,
    30.0,
    35.0,
    40.0,
]
DIFF_LEVELS = [
    -15.0,
    -10.0,
    -5.0,
    -2.0,
    -1.0,
    -0.5,
    -0.2,
    0.2,
    0.5,
    1.0,
    2.0,
    5.0,
    10.0,
    15.0,
]


def main() -> int:
    param = CoreParameter()
    param.no_viewer = True
    param.reference_data_path = REFERENCE_DATA_PATH
    param.test_data_path = TEST_DATA_PATH
    param.results_dir = RESULTS_DIR
    param.case_id = "CRU_IPCC"
    param.sets = ["lat_lon"]
    param.variables = ["TREFHT"]
    param.seasons = ["ANN"]
    param.regions = ["land_60S90N"]
    param.regrid_method = "bilinear"
    param.multiprocessing = False
    param.num_workers = 8
    param.main_title = "TREFHT ANN land_60S90N"
    param.contour_levels = CONTOUR_LEVELS
    param.test_name = "v3.LR.historical_0051"
    param.short_test_name = "v3.LR.historical_0051"
    param.ref_name = "CRU"
    param.reference_name = "CRU Global Monthly Mean T Land"
    param.short_ref_name = "CRU"
    param.diff_title = "Model - Observations"
    param.diff_levels = DIFF_LEVELS

    runner.sets_to_run = ["lat_lon"]
    logger.info("Running lat_lon QA case for issue 1044")
    runner.run_diags([param])

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
