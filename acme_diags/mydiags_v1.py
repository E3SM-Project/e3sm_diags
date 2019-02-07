reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'

test_data_path = '/Users/zhang40/Documents/ACME_simulations/'

test_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'

backend = 'cartopy'
diff_title = 'Model - Obs.'

results_dir = 'table_all'
seasons = ['ANN','JJA']

sets = ['lat_lon']

#sets = ['5', '13']

#def albedo_obs(rsdt, rsut):
#    """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
#    var = rsut / rsdt
#    var.units = "dimensionless"
#    var.long_name = "TOA albedo"
#    return var
#
#derived_variables = {
#    'ALBEDO': {
#        ('rsdt', 'rsut'): albedo_obs
#    }
#}
