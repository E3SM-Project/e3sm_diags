import os
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.run import runner

simulations = ['tc-v2.LR.amip.2000_2014']
sim_names = ['v2.LR.amip_0101']

data_path = '/global/homes/c/chengzhu/tests/tc_analysis/'

for idx, sim in enumerate(simulations):
    print(sim)
    param = CoreParameter()
    param.test_data_path = data_path+sim
    param.test_name = sim_names[idx]
    param.test_start_yr = "2000"
    param.test_end_yr = "2014"

    param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/tc-analysis'
    param.ref_start_yr = '1979'
    param.ref_end_yr = '2018'

    prefix = f'/global/cfs/cdirs/e3sm/www/chengzhu/tc_analysis_test/'
    param.results_dir = os.path.join(prefix, sim)
    runner.sets_to_run = ['tc_analysis']
    runner.run_diags([param])
