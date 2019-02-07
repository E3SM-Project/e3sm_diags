reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
#reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
test_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

#test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
#test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01'
#test_name = 'AMWG obs'

backend = 'cartopy'
diff_title = 'New obs. - Old Obs.'
test_title = 'New obs'
reference_title = 'Old obs'
diff_colormap = 'bl_to_darkred'
results_dir = 'newobs_sst_ceres_erai_merra2_vs_oldobs'

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
