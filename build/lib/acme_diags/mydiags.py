reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
test_name = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01'

backend = 'vcs'
diff_title = 'Model - Obs.'
diff_colormap = 'bl_to_darkred'
results_dir = 'acme_diags_demo_vcs210'

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
