import os

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

param = CoreParameter()

#param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology'
#param.test_data_path = '/global/cfs/cdirs/e3sm/zhang40/e3sm_diags_for_EAMxx/data/Cess'
param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology'
param.test_data_path = '/global/cfs/cdirs/e3sm/chengzhu/eamxx/post/data/rgr'
param.test_name = 'eamxx_decadal'
param.seasons = ["ANN"]
#param.save_netcdf = True

prefix = '/global/cfs/cdirs/e3sm/www/zhang40/tests/eamxx'
param.results_dir = os.path.join(prefix, 'eamxx_decadal_1212')

runner.sets_to_run = ["lat_lon",
        "zonal_mean_xy",
        "zonal_mean_2d",
        "zonal_mean_2d_stratosphere",
        "polar",
        "cosp_histogram",
        "meridional_mean_2d",
        "annual_cycle_zonal_mean",
        ]

runner.run_diags([param])

