import e3sm_diags
from e3sm_diags.driver import utils
import cdms2
import cdutil
import numpy

import csv

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
#mwso4 = 115.0
#mws = 32.066
#units_conv_factor_S = units_conv_factor/mwso4*mws # kg/s to TgS/yr


species = ["bc", "dst", "mom", "ncl","pom","so4","soa"]
SPECIES_NAMES = {"bc": "Black Carbon",
    "dst": "Dust",
    "mom": "Marine Organic Matter",
    "ncl": "Sea Salt",
    "pom": "Primary Organic Matter",
    "so4": "Sulfate",
    "soa": "Secondary Organic Aerosol"}
missing_value = 999.999
metrics_dict = {}

seasons = ["ANN"]

for season in seasons:
    for aerosol in species:
        print(f'Aerosol species: {aerosol}')
        wetdep = test_data.get_climo_variable(f'{aerosol}_SFWET', season) 
        drydep = test_data.get_climo_variable(f'{aerosol}_DDF', season) 
        srfemis = test_data.get_climo_variable(f'SF{aerosol}', season) 
        mass = test_data.get_climo_variable(f'Mass_{aerosol}', season) 
        area = test_data.get_extra_variables_only(
                    f'{aerosol}_DDF', season, extra_vars=["area"]
                )
        hyam, hybm, ps = test_data.get_extra_variables_only(
            f'Mass_{aerosol}', season, extra_vars=["hyai", "hybi", "PS"]
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
        burden_total= numpy.sum(numpy.sum(burden*area_m2, axis = 0), axis=0)*1e-9 # kg to Tg
        print(f'{aerosol} Burden (Tg): ',f'{burden_total:.3f}')
        sink = numpy.sum(numpy.sum((drydep-wetdep)*area_m2,axis = 0), axis=0)*units_conv_factor
        drydep = numpy.sum(numpy.sum((drydep)*area_m2,axis = 0), axis=0)*units_conv_factor
        wetdep = numpy.sum(numpy.sum((-wetdep)*area_m2,axis = 0), axis=0)*units_conv_factor
        srfemis = numpy.sum(numpy.sum((srfemis)*area_m2,axis = 0), axis=0)*units_conv_factor
        print(f'{aerosol} Sink (Tg/year): ',f'{sink:.3f}')
        print(f'{aerosol} Lifetime (days): ',f'{burden_total/sink*365:.3f}')
        metrics_dict[aerosol] = {
    	"Surface Emission (Tg/yr)": f'{srfemis:.3f}',
    	"Sink (Tg/s)": f'{sink:.3f}',
    	"Dry Deposition (Tg/yr)": f'{drydep:.3f}',
    	"Wet Deposition (Tg/yr)": f'{wetdep:.3f}',
    	"Burden (Tg)": f'{burden_total:.3f}',
    	"Lifetime (Days)": f'{burden_total/sink*365:.3f}',
        }
        with open(f'aerosol_table_{season}.csv', "w") as table_csv:
            writer = csv.writer(
                table_csv,
                delimiter=",",
                quotechar="'",
                quoting=csv.QUOTE_MINIMAL,
                lineterminator='\n',
            )
            #writer.writerow([" ", "test","ref",])
            for key, values in metrics_dict.items():
                writer.writerow([SPECIES_NAMES[key]])
                print('key',key, values)
                for value in values:
                    print(value)
                    line = []
                    line.append(value)
                    line.append(values[value])
                    writer.writerows([line])
                writer.writerows([""])




