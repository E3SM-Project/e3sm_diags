from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
from e3sm_diags.run import runner
import sys

param = ZonalMean2dParameter()

param.reference_data_path = '/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/climatology'
param.test_data_path = '/lcrc/group/e3sm/ac.zhang40/example_v3/v3.LR.historical_0051/post/atm/180x360_aave/clim/30yr'
param.results_dir = '/lcrc/group/e3sm/public_html/diagnostic_output/ac.tvo/tests/pdf_size_1'
param.case_id = 'ERA5'
param.run_type = 'model_vs_obs'
param.sets = ['zonal_mean_2d']
param.variables = ['T']
param.seasons = ['ANN']
param.regions = ['global']
param.plevs = [50.0, 100.0, 150.0, 200.0, 250.0, 300.0, 350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, 750.0, 800.0, 850.0, 900.0, 950.0, 1000.0]
param.multiprocessing = True
param.num_workers = 8
param.main_title = 'T ANN global'
param.backend = 'cartopy'
param.output_format = ['png']
param.canvas_size_w = 1212
param.canvas_size_h = 1628
param.figsize = [8.5, 11.0]
param.dpi = 150
param.arrows = True
param.contour_levels = [180, 185, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 295, 300]
param.test_name = 'v3.LR.historical_0051'
param.short_test_name = 'v3.LR.historical_0201'
param.test_colormap = 'cet_rainbow.rgb'
param.ref_name = 'ERA5'
param.reference_name = 'ERA5 Reanalysis'
param.reference_colormap = 'cet_rainbow.rgb'
param.diff_title = 'Model - Observations'
param.diff_colormap = 'diverging_bwr.rgb'
param.diff_levels = [-7, -6, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 6, 7]
param.granulate = ['variables', 'seasons', 'regions']
param.selectors = ['sets', 'seasons']
param.save_netcdf = True
param.output_format_subplot = ['pdf']

# cfg_filepath = 'auxiliary_tools/debug/987-pdf-no-plots/mvce.cfg'
# sys.argv.extend(['-d', cfg_filepath])
runner.sets_to_run = ['zonal_mean_2d']
runner.run_diags([param])