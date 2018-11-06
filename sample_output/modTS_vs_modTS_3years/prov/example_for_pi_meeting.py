short_test_name = 'historical_H1'
short_ref_name = 'historical_H1'

test_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'
reference_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'
run_type = 'model_vs_model'

test_timeseries_input = True
test_start_yr = '2011'
test_end_yr = '2013'

ref_timeseries_input = True
ref_start_yr = '1850'
ref_end_yr = '1852'

#sets = ['lat_lon']
#seasons=['ANN']

results_dir = '/var/www/acme/acme-diags/shaheen2/dataset/modTS_vs_modTS_3years'
diff_title = 'Model (2011-2013) - Model (1850-1852)'

multiprocessing = True
num_workers = 64

