#This script is to combine precipication from GPCP and evaporation from OAFlux to derive P minus E variable
#Generated Oct25th by Jill Zhang 

import cdms2 
import os
import numpy

# Set nc classic as outputs
cdms2.setCompressionWarnings(0) ; # Suppress warnings
cdms2.setNetcdfShuffleFlag(0)
cdms2.setNetcdfDeflateFlag(1) ; # was 0 130717
cdms2.setNetcdfDeflateLevelFlag(9) ; # was 0 130717
cdms2.setAutoBounds(1) ; # Ensure bounds on time and depth axes are generated

input_path1="/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/GPCP/time_series/"               #location of input data
input_path2="/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/OAFlux/time_series/"               #location of input data
filename1='PRECT_197901_201712.nc'
filename2='QFLX_197901_201312.nc'

fin1= cdms2.open(input_path1+filename1) 
fin2= cdms2.open(input_path2+filename2) 

startyear='1979'  #start year of the dataset
endyear='2013'    #end year of the dataset

prect = fin1('PRECT',time=("1979-01-15","2013-12-15", 'ccb'))
qflux = fin2('QFLX',time=("1979-01-15","2013-12-15", 'ccb'))
#print(prect.getTime().asComponentTime())
#print(qflux.getTime().asComponentTime())

qflux_reg = qflux.regrid(prect.getGrid(), regridTool='esmf', regridMethod='linear')

PminusE=prect/24.0/3600.0 - qflux_reg #conver to kg/m2/s
PminusE.setAxis(0,prect.getTime())
PminusE.id = 'PminusE'
PminusE.units = 'kg/m2/s'
PminusE.long_name = 'Precipitation - evaporation rate'
PminusE.standard_name = 'P - E flux'
PminusE._FillValue = 1.0e+20
PminusE.missing_value = 1.0e+20

delattr(PminusE,'least_significant_digit')
delattr(PminusE,'level_desc')
delattr(PminusE,'actual_range')
delattr(PminusE,'precision')
delattr(PminusE,'parent_stat')
delattr(PminusE,'dataset')
delattr(PminusE,'var_desc')
delattr(PminusE,'valid_range')

output_path = '/p/user_pub/e3sm/zhang40/analysis_data_e3sm_diags/GPCP_OAFLux/time_series/'
output_filename='PminusE_197901_201312.nc'
fout = cdms2.open(output_path+output_filename,'w')
fout.write(PminusE)
year_str="".join([startyear,'-',endyear])
setattr(fout,'yrs_averaged',year_str)
setattr(fout,'data_name','P - E computed using precipitation rate from GPCP (1979-2013) and evaporation rate from OAFlux (1979-2013)')
print(fout.attributes)
fout.close()


