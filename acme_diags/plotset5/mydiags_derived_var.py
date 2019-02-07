reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
test_name = '20160520.A_WCYCL1850.ne30'

backend = 'vcs'
diff_title = 'Test - Reference'
diff_colormap = 'bl_to_darkred'

def albedo_obs(rsdt, rsut):
    """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
    var = rsut / rsdt
    var.units = "dimensionless"
    var.long_name = "TOA albedo"
    return var

derived_variables = {
    'ALBEDO': {
        ('rsdt', 'rsut'): albedo_obs
    }
}
