import e3sm_diags
from e3sm_diags.driver import utils
import cdms2
import cdutil
import numpy

from e3sm_diags.parameter.core_parameter import CoreParameter

param = CoreParameter()
param.test_name = 'v2.LR.historical_0101'
#param.test_name = 'v2.LR.F2010'
param.test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
test_data = utils.dataset.Dataset(param, test=True)

rearth = 6.37122e6 #km
units_conv_factor = 86400.0*365.0*1e-9 # kg/s to Tg/yr
# TODO:
# Convert so4 unit to TgS

species = ["bc", "dst", "mom", "ncl","pom","so4","soa"]
species_long_name = ["Black Carbon", "Dust", "Marine Organic Matter", "Sea Salt", "Primary Organic Matter","Sulfate", "Secondary Organic Aerosol"]

for aerosol in species:
    print(f'Aerosol species: {aerosol}')
    wetdep = test_data.get_climo_variable(f'{aerosol}_SFWET', 'ANN') 
    drydep = test_data.get_climo_variable(f'{aerosol}_DDF', 'ANN') 
    srfemis = test_data.get_climo_variable(f'SF{aerosol}', 'ANN') 
    mass = test_data.get_climo_variable(f'Mass_{aerosol}', 'ANN') 
    area = test_data.get_extra_variables_only(
                f'{aerosol}_DDF', 'ANN', extra_vars=["area"]
            )
    hyam, hybm, ps = test_data.get_extra_variables_only(
        f'Mass_{aerosol}', 'ANN', extra_vars=["hyai", "hybi", "PS"]
    )
    area_m2 = area * rearth**2
    
    
    p0 = 100000.0#1000.0  # mb
    ps = ps #/ 100.0  # convert unit from 'Pa' to mb
    pressure_levs = cdutil.vertical.reconstructPressureFromHybrid(ps, hyam, hybm, p0)
    pressure_levs.units = "Pa"
    
    #(72,lat,lon)
    delta_p = numpy.diff(pressure_levs,axis = 0)
    mass_3d = mass*delta_p/9.8 #mass density * mass air   kg/m2
    #print(aero_3d.shape)
    burden = numpy.nansum(mass_3d,axis = 0)   #kg/m2
    burden_total= numpy.sum(numpy.sum(burden*area_m2, axis = 0), axis=0)*1e-9
    print(f'{aerosol} Burden (Tg): ',f'{burden_total:.3f}')
    sink = numpy.sum(numpy.sum((drydep-wetdep)*area_m2,axis = 0), axis=0)*units_conv_factor
    print(f'{aerosol} Sink (Tg/year): ',f'{sink:.3f}')
    print(f'{aerosol} Lifetime (days): ',f'{burden_total/sink*365:.3f}')

