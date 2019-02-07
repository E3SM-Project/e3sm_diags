import json

case_id = 'set5_ANN_PRECT_GPCP'
#reference_data_path = '/space1/test_data/obs_for_diagnostics/'
#test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'
#reference_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'
#test_data_path = '/Users/shaheen2/github/acme_diags/acme_diags/'

reference_data_path = '/Users/zhang40/Documents/obs_for_diagnostics/'
test_data_set = '20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model
test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
#reference_data_set = 'GPCP_ANN_climo.nc'  # observation

## List of variables in NCAR name convention:
#variables = ['PRECT','TREFHT','SST','PRECT','PREH2O','QFLX','LHFLX','SHFLX','TGCLDLWP','FLNS','FSNS','ALBEDO','ALBEDOC','FLUT','FLUTC','FSNTOA','FSNTOAC','SWCF','LWCF','FLDS','FLDSC','FSDS','FSDSC','FLNS','FLNSC','FSNS','FSNSC','LWCFSRF','SWCFSRF','CLDHGH','CLDLOW','CLDMED','CLDTOT','MEANPTOP','MEANTTOP','MEANTAU','TCLDAREA','PSL','T','U','Z']

# read in json file for observationa information
filename='/Users/zhang40/Documents/AIMS/repo/acme-diags/acme_diags/plotset5/obs_info_dictionary.json'

with open(filename) as json_data:
   obs_data  = json.load(json_data)



# variables and their observational data name
variables = {'PRECT': ['ERAI','GPCP','LEGATES'],
              'T': ['EARI','MERRA']}

## Regions for metrics calculation (var specific)
regions = {'PRECT' : [None,'land','ocean', 'TROPICS'],
           'TREFHT' : [None,'land','ocean'],
           'TGCLDLWP' : [None,'ocean']} 

## Pressure level slice for 3D variable 
plev = {'T' : [850.,200.],
        'U' : [200.],
        'Z' : [300.,200.]} 

## Seasons
season =['ANN','DJF','MAM','JJA','SON']
#season =['ANN','DJF','MAM','JJA','SON','01','02','03','04','05','06','07','08','09','10','11','12']

## Regrid options
regrid_tool = 'esmf'
regrid_method = 'linear'

## Observational data dictionary:
set5_observations = 'set5_obs_info_dictionary.json'


## Output file
ext = '.png' 

if var == ['T','U','Z']:
    main_title = var + plev + season
  output_file = "%(variable)_%(plev)_%(season)_%(regridMethod)_%(regions)_%(ext)"

else:
    main_title = var + season
  output_file = "%(variable)_%(season)_%(regridMethod)_%(regions)_%(ext)"


## More plotting setting
#test_name = '1850_alpha6_01 (yrs0070-0099)'
#test_title = 'Model'
#test_colormap = ''
##test_levels = [0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
#test_levels = []
#test_units = 'mm/day'
#
#reference_name = 'GPCP (yrs1979-2009)'
#reference_title = 'Observation'
#reference_colormap = ''
##reference_levels = [0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
#reference_levels = []
#reference_units = 'mm/day'
#
#diff_name = ''
#diff_title = 'Model - Observation'
#diff_colormap = 'bl_to_darkred'
##diff_levels = [-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6]
#diff_levels = []
#diff_units = 'mm/day'

canvas_size_w = 1212
canvas_size_h = 1628
