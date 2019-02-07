reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
#test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01'
test_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'

backend = 'cartopy'
diff_title = 'Model - Obs.'
diff_colormap = 'rainbow'
test_colormap = 'rainbow'
reference_colormap = 'bl_to_darkred'
results_dir = 'tier1b_cloud'
contour_levels =[0,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.]
diff_levels = [-20,-15,-10,-5,0,5,10,15,20]
#seasons = ['ANN','DJF','MAM','JJA','SON']
seasons = ['ANN']

sets = [5]

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
