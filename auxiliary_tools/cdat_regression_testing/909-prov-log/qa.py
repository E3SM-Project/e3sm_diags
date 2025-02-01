import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.run import Run, runner

simulations = [
    "tc-v1.HR.0026_0035",
    "tc-v2.LR.2000_2014",
    "tc-v3.HR.0006_0025",
    "tc-v3.LR.2000_2014",
]
sim_names = [
    "theta.20180906.branch_noCNT.A_WCYCL1950S_CMIP6_HR.ne120_oRRS18v3_ICG",
    "v2.LR.historical_0101",
    "20240609.piCtl.ne120pg2_r025_RRSwISC6to18E3r5.chrysalis.test1",
    "extendedOutput.v3.LR.historical_0101",
]

data_path = "/global/homes/c/chengzhu/tests/tc_analysis/"

for idx, sim in enumerate(simulations):
    print(sim)
    # runner = Run()

    param = CoreParameter()
    param.multiprocessing = True
    param.test_data_path = data_path + sim
    param.test_name = sim_names[idx]
    param.test_start_yr = sim.split(".")[-1][0:4]
    param.test_end_yr = sim.split(".")[-1][5:9]

    param.reference_data_path = (
        "/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/tc-analysis"
    )
    param.ref_start_yr = "1979"
    param.ref_end_yr = "2018"

    prefix = "/global/cfs/cdirs/e3sm/www/vo13/tc_analysis_test/"
    param.results_dir = os.path.join(prefix, sim)
    runner.sets_to_run = ["tc_analysis"]
    runner.run_diags([param])
