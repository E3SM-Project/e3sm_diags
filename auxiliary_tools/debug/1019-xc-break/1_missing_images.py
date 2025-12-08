from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner

# Create a CoreParameter object and set its attributes
parameters = CoreParameter()
parameters.no_viewer = True
parameters.reference_data_path = "/lcrc/group/e3sm/diagnostics/observations/Atm/climatology"
parameters.test_data_path = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis/climatology/rgr/"
parameters.results_dir = "/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/annual_cycle_xr_2025.11.0"
parameters.case_id = "CERES-EBAF-surface-v4.1"
parameters.sets = ["annual_cycle_zonal_mean"]
parameters.variables = ["FSNSC"]
parameters.multiprocessing = True
parameters.num_workers = 24
parameters.main_title = "FSNSC ANNUALCYCLE global"
parameters.contour_levels = [0, 50, 100, 150, 200, 250, 300, 350, 400]
parameters.test_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
parameters.ref_name = "ceres_ebaf_surface_v4.1"
parameters.reference_name = "CERES-EBAF v4.1"
parameters.diff_levels = [-75, -50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50, 75]

# Pass the CoreParameter object to run_diags
runner.sets_to_run = ['annual_cycle_zonal_mean']
runner.run_diags([parameters])