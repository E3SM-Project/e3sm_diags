# Location of the data.
reference_data_path = '/Users/zhang40/Documents/ACME_simulations/obs_climo/climatology/'
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
# Name of the test model data, used to find the climo files.
test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6'
# An optional, shorter name to be used instead of the test_name.
short_test_name = 'WCYCL1850.ne30'

# What plotsets to run the diags on.
sets = ['meridional_mean_2d']
seasons = ['ANN']
# Name of the folder where the results are stored.
results_dir = 'model_to_obs_meridional_mean'

# Below are more optional arguments.

# 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
backend = 'mpl'
# Title of the difference plots.
diff_title = 'Model - Obs.'

