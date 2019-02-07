#!/usr/bin/env python
import numpy
import cdutil
import cdms2
import acme_diags.acme_parser
import plot_set_5

    

parser = acme_diags.acme_parser.ACMEParser()
parameter = parser.get_parameter()

reference_data_path = parameter.reference_data_path
reference_data_set = parameter.reference_data_set # observation
test_data_path = parameter.test_data_path
test_data_set = parameter.test_data_set # model

var = parameter.variables
season = parameter.season


f_obs = cdms2.open(reference_data_path + reference_data_set)
f_mod = cdms2.open(test_data_path + test_data_set)

if var == 'PRECT':
    obs_pr = f_obs(var, longitude=(-180, 540))
    mod_pr = (f_mod('PRECC', longitude=(-180, 540)) + f_mod('PRECL', longitude=(-180, 540)))*3600.0*24.0*1000.0
    mod_pr.units = 'mm/day'
    obs_name = 'GPCP (yrs1979-2009)'
elif var == 'T':
    obs_pr = f_obs(var)#, longitude=(-180, 180))
    mod_pr = f_mod(var)#, longitude=(-180, 180))
    print 'aaa',mod_pr[:,:,0,0]
    #obs_name = 'ECMWF (yrs unknown)'
    # set the level to compare, this should be a parameter
    plev = 220 #mb

if mod_pr.ndim == 4: # var(time,lev,lon,lat) 
    obs_plv = obs_pr.getLevel()
    mod_plv = mod_pr.getLevel()
    print mod_plv
    
    if mod_plv.long_name.lower().find('hybrid') != -1: # var(time,lev,lon,lat) convert from hybrid level to pressure
        hyam = f_mod('hyam')
        hybm = f_mod('hybm')
        ps = f_mod('PS')/100.    #convert unit from 'Pa' to mb
        p0 = 1000. #mb
        levels_orig = cdutil.vertical.reconstructPressureFromHybrid(ps,hyam,hybm,p0)
        levels_orig.units = 'mb'

        #if plev in levels_orig[:].tolist():
        print 'lll',levels_orig[:,:,0,0]
        mod_pr_p=cdutil.vertical.logLinearInterpolation(mod_pr, levels_orig, plev)
        print 'levels_orig', levels_orig.shape 
        print levels_orig[0,:,0,0],mod_pr[0,:,0,0]
        #plv17 = [10.,30.,50.,70.,100.,150.,220.,250.,300.,400.,500.,
        #                 600.,700.,775.,850.,925.,1000.]
        #mod_pr_p=cdutil.vertical.logLinearInterpolation(mod_pr, levels_orig, plv17)
        #mod_plv_p=mod_pr_p.getLevel()
        #print mod_plv_p[:]
        #levels_orig_p=mod_pr_p.clone()
        #for ilev in range(len(mod_plv_p[:])):
        #    levels_orig_p.data[:,ilev,:,:]=float(mod_plv_p[ilev])
        #mod_pr_p_p=cdutil.vertical.logLinearInterpolation(mod_pr_p, levels_orig_p, plev)
        #print 'ccc',mod_pr_p_p[:,:,0,0],mod_pr_p[:,:,0,0]
        ##else:
        print mod_pr_p.shape

    else: # mod_plv.long_name.lower().find('pressure') !=-1: 
        levels_orig = mod_pr[:,:,:,:].clone()#.clone()
        levels_orig.units = 'mb'
        levels_orig.id = 'pressure'
        print levels_orig.shape
        for ilev in range(len(mod_plv[:])):
            #levels_orig[:,ilev,:,:]=mod_plv[len(mod_plv[:])-1-ilev]
            levels_orig.data[:,ilev,:,:]=mod_plv[ilev]
        mod_pr_p=cdutil.vertical.logLinearInterpolation(mod_pr[:,::-1,:,:], levels_orig[:,::-1,:,:], plev)
        print 'levels_orig', levels_orig.shape 
        print levels_orig[0,:,0,0],mod_pr[0,:,0,0]
        
        print 'ccc',mod_pr_p[:,:,0,0],mod_pr[:,:,0,0]
        
#        quit()
        #mod_pr_p = mod_pr.pressureRegrid([225])
        #mod_pr_p = mod_pr

#    else: 
#        print( 'Vertical level is neither hybrid nor pressure.')
#        quit()
        
    # Create pressure grid to interpolate reference data
    # Save this for set4 metrics
    #obs_levels_orig = cdms2.createAxis(plv17,id = 'lev')
    #obs_pr_p = obs_pr.pressureRegrid(obs_levels_orig)
    #obs_pr = obs_pr_p[:,plev_ind,:,:]
    try: 
        #plev_ind = plv17.index(plev)
        #mod_pr = mod_pr_p[:,plev_ind,:,:]
        mod_pr = mod_pr_p#[:,plev_ind,:,:]
        print mod_pr.shape
    except:
        print( 'Spicified level is not a reference level for model')
    try:
        plev = 200
        plev_ind = obs_plv[:].tolist().index(plev)
        #obs_pr = obs_pr[:,plev_ind:plev_ind+1,:,:]
        obs_pr = obs_pr[:,plev_ind,:,:]
        print obs_pr.shape
    except:
        print( 'Spicified level is not a reference level')

axes1 = mod_pr.getAxisList()
axes2 = obs_pr.getAxisList()

# For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
if len(axes1[1]) <= len(axes2[1]): # use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
    model_grid = mod_pr.getGrid()
    mod_pr_reg = mod_pr
    obs_pr_reg = obs_pr.regrid(model_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)
else:
    obs_grid = obs_pr.getGrid()
    obs_pr_reg = obs_pr
    mod_pr_reg = mod_pr.regrid(obs_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)

if var == 'T' and plev == 850:
    levels = [230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300]
    diff_levels = [-8, -6, -5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5, 6, 8]
    parameter.reference_levels = levels
    parameter.test_levels = levels
    parameter.diff_levels = diff_levels
elif var == 'T' and plev == 200:
    levels = [190, 193, 196, 199, 202, 205, 208, 211, 214, 217, 220, 223, 226, 229, 232]
    diff_levels = [-10, -8, -6, -4, -3, -2, -1, 0, 1, 2, 3, 4, 6, 8, 10]
    parameter.reference_levels = levels
    parameter.test_levels = levels
    parameter.diff_levels = diff_levels

if var == 'T':
    parameter.main_title = ' '.join([var, str(plev), 'mb', season])

plot_set_5.plot(obs_pr, mod_pr, obs_pr_reg, mod_pr_reg, parameter)
