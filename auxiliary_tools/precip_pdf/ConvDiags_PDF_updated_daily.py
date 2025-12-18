#!/usr/local/uvcdat/1.3.1/bin/python
"""
Created to produce netcdfs of daily-mean PRECT and PRECC on a 1degx1deg 
grid, consistent with the GPCP1DD dataset. 

Input:
As inputs, the model takes a netcdf that contains latxlon gridded PRECT (PRECC)
output from the model. 

Output:
The output will be a netcdf that includes two 3D variables (lat,lon,precip rate) containing 
the frequency PDF (FREQPDF) and amount PDF (AMNTPDF) for each precipitation rate bin 
at each 1degx1deg box. bounds_binvalue (edges of each precip bin). bincenter gives the arithmetic
mean of the bounds values for each bin. 

Global/regional mean PDFs can be reconstructed from the gridded PDFs (see simple plotting below). 

@author: terai1
Original date: 24 April 2015 (adapted from original tier1b diagnostics)
Modified date: 2 Nov 2020
"""
import argparse,datetime,gc,re,sys,time
import cdms2 as cdm
import MV2 as MV     #stuff for dealing with masked values.
import cdtime        #for obtaining time values for splicing data
import glob
import os
import cdutil
import genutil   #for averager
import numpy as np
import WC_diag_amwg

cdm.setAutoBounds('on')  #in case timebounds is not set


#  **************************************************************
#     
#    USERS to modify and check
#      1 Specify location of model data (daily mean regridded data)
#      2 Specify location for figure output
#      3 Specify which precipitation variable to create PDFs of (PRECT,PRECC,PRECL)
#  
#  **************************************************************


model_hostpath='/global/cfs/cdirs/e3sm/user/location_ofregriddedh3data'   #'/directory/regridded_model_data/'
model_hostpath='/global/cfs/cdirs/e3sm/chengzhu/tests/zppy_example_v3/v3.LR.amip_0101/post/atm/180x360_aave/ts/daily/10yr/'
obs_hostpath='/global/cfs/cdirs/e3sm/terai/Obs_datasets/GPCP_PDF/'          # all users associated with e3sm group should be able to access
figure_directory='./Figures/'
variable_2_runwith='PRECT'


PrecipPDF_YorN=1    #Yes is 1; No is 0

#   **************************************************************

print(model_hostpath)





if PrecipPDF_YorN>0:
    #Access GPCP netcdf (mainly to obtain grid)
    
    f_in=cdm.open(''.join([obs_hostpath,'GPCP_1DD_v1.2_199610_201312.nc']))
    obs_freq_pdf=f_in['PREC']
    
    #   **********************************************************************
    #   identify year for which to take data from - here I've specified 1981 for comparison
    #   **********************************************************************
    fis=glob.glob(os.path.join(model_hostpath,''.join(['PRECT_*.nc']))) #find gridded data - match PRECT files
    
    data=[] #initialize lists to put stuff in
    locs=[]
    cnt=0
    print('Got past fis')
    print(fis)
    for fi in fis:  #Create a netcdf with the pdfs for each model netcdf file
        variablename=variable_2_runwith
        f_model_dailyprecip=cdm.open(fi)
        print('Loaded model data')
        model_handle=f_model_dailyprecip
        case= model_handle.attributes['case']
        
        pdfpath=model_hostpath
        model_filename_str=fi
        # Parse filename: PRECT_199501_200412.nc
        basename = os.path.basename(fi)
        name_parts = basename.replace('.nc', '').split('_')
        model_time_str = '_'.join(name_parts[1:])  # Get "199501_200412"
        model_pdf_filename='_'.join([case,model_time_str,variablename,'PDF.nc'])
        filename_check=os.path.join(pdfpath,model_pdf_filename)
        
        
        if os.path.isfile(filename_check):
            print("netcdf file containing PDF already exists in path")
        else:
            mv=f_model_dailyprecip(variablename)
            # Skip regridding - data already on 180x360_aave grid
            #gpcp_grid=obs_freq_pdf.getGrid()
            #mv=mv.regrid(gpcp_grid,regridTool='esmf')#,regridMethod='conserve')
            mv.id=variablename
            #If the units are not in mm/day, convert units
            target_units_prect='mm/day'
            if mv.units=='kg m-2 s-1' or mv.units=='kg/m2/s' or mv.units=='kg m^-2 s^-1':
                mv.units='mm/s'
            ## Change units of precip from m s-1 to mm/day
            mv=mv*3600.*24.*1000.
            #mv=mv*24. #change from mm/hr to mm/day
            mv.units=target_units_prect
            mv.id=variablename
            WC_diag_amwg.create_precip_PDF_netcdf(mv, model_handle, model_pdf_filename, path=pdfpath,model_time=model_time_str)
            print("---------- pdf created")




