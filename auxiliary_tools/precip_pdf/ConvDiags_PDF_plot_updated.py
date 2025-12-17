#!/usr/local/uvcdat/1.3.1/bin/python
"""
Plots PDFs of precipitation rate

@author: terai1
"""
#Import necessary libraries

import cdms2
import cdutil
import genutil
import matplotlib.pyplot as plt
import numpy as np
import requests
import WC_diag_amwg
cdms2.setAutoBounds(1)


# Access PDF files

# ***** NOTE: The script below is for 3-hourly data, you can find the daily PDFs here:
# Model: /global/cscratch1/sd/terai/E3SM_NGD-ATMPhys/E3SM_alpha5_59/wMG2/daily_rgr/
# Obs: /global/cfs/cdirs/e3sm/terai/Obs_datasets/TRMM_PDF/TRMM_3B42_GPCPgrid_1998-2007_PDF.nc
# GPCP: /global/cfs/cdirs/e3sm/terai/Obs_datasets/GPCPv1pt2_PREC_pdf.nc
# **********************************

# Model PDF files
data_directory='/global/cfs/cdirs/e3sm/user/location_ofregriddedh3data'   #user modify
data_filename='CASENAME_2011_PRECT_PDF.nc'   #user modify
f_in_MG2_2011=cdms2.open(''.join([data_directory,data_filename]))
data_filename='CASENAME_2012_PRECT_PDF.nc'   #user modify
f_in_MG2_2012=cdms2.open(''.join([data_directory,data_filename]))
data_filename='CASENAME_2013_PRECT_PDF.nc'   #user modify
f_in_MG2_2013=cdms2.open(''.join([data_directory,data_filename]))



#TRMM PDF files
TRMM_PDF_file='/global/cfs/cdirs/e3sm/terai/Obs_datasets/TRMM_PDF/Three_hourly/TRMM.3B42v7_2000_PRECT_PDF.nc'
f_in_TRMM_PDF_2000=cdms2.open(TRMM_PDF_file)
TRMM_PDF_file='/global/cfs/cdirs/e3sm/terai/Obs_datasets/TRMM_PDF/Three_hourly/TRMM.3B42v7_2001_PRECT_PDF.nc'
f_in_TRMM_PDF_2001=cdms2.open(TRMM_PDF_file)
TRMM_PDF_file='/global/cfs/cdirs/e3sm/terai/Obs_datasets/TRMM_PDF/Three_hourly/TRMM.3B42v7_2002_PRECT_PDF.nc'
f_in_TRMM_PDF_2002=cdms2.open(TRMM_PDF_file)


# Retrieve Tropical mean PDF of FREQ and AMNT PDFs
bincenter,MG2_2011_trop_PRECTFREQPDF,MG2_2011_trop_PRECTAMNTPDF=WC_diag_amwg.Tropical_PDF_values(f_in_MG2_2011,'PRECT')
bincenter,MG2_2012_trop_PRECTFREQPDF,MG2_2012_trop_PRECTAMNTPDF=WC_diag_amwg.Tropical_PDF_values(f_in_MG2_2012,'PRECT')
bincenter,MG2_2013_trop_PRECTFREQPDF,MG2_2013_trop_PRECTAMNTPDF=WC_diag_amwg.Tropical_PDF_values(f_in_MG2_2013,'PRECT')


bincenter,TRMM_2000_trop_PRECTFREQPDF,TRMM_2000_trop_PRECTAMNTPDF  =WC_diag_amwg.Tropical_PDF_values(f_in_TRMM_PDF_2000,'PRECT')
bincenter,TRMM_2001_trop_PRECTFREQPDF,TRMM_2001_trop_PRECTAMNTPDF  =WC_diag_amwg.Tropical_PDF_values(f_in_TRMM_PDF_2001,'PRECT')
bincenter,TRMM_2002_trop_PRECTFREQPDF,TRMM_2002_trop_PRECTAMNTPDF  =WC_diag_amwg.Tropical_PDF_values(f_in_TRMM_PDF_2002,'PRECT')


MG2_array=np.zeros((129,3))
MG2_array[:,0]=MG2_2011_trop_PRECTFREQPDF
MG2_array[:,1]=MG2_2012_trop_PRECTFREQPDF
MG2_array[:,2]=MG2_2013_trop_PRECTFREQPDF


TRMM_array=np.zeros((129,3))
TRMM_array[:,0]=TRMM_2000_trop_PRECTFREQPDF
TRMM_array[:,1]=TRMM_2001_trop_PRECTFREQPDF
TRMM_array[:,2]=TRMM_2002_trop_PRECTFREQPDF

#Plot the model and obs global mean frequency pdf'                                                                              
fig = plt.figure(figsize=(9,8))
for i in np.arange(3):
    plt.plot(bincenter,MG2_array[:,i],color='tab:blue')
    plt.plot(bincenter,TRMM_array[:,i],color='tab:brown')
    
plt.plot(bincenter,np.mean(MG2_array,axis=1),linewidth=2,color='tab:blue',label='CASENAME')    #user modify

plt.plot(bincenter,np.mean(TRMM_array,axis=1),color='tab:brown',linewidth=2,label='TRMM')
plt.ylabel('df/dlog(P)',fontsize=16)
plt.xlabel('P (mm/day)',fontsize=16)
plt.xscale('log')
plt.yticks(size=16)
plt.xticks(size=16)
plt.legend(loc='upper right',fontsize=16)
plt.title('Tropical-mean (30$^{\circ}$S-30$^{\circ}$N) Precipitation Frequency PDF (3-hourly)',fontsize=16)
plt.savefig('Figurename')   #user modify

# for extreme precip
bounds_binvalue=f_in_MG2_2011('bounds_binvalue',binvalue=[0.1,600.])
bounds_log10diff=np.log10(np.array(bounds_binvalue)[0,1]/np.array(bounds_binvalue)[0,0])
bounds_lindiff=np.array(bounds_binvalue)[:,1]-np.array(bounds_binvalue)[0,0] # To scale the frequency diagram

#Average across three years and also scale to allow for linear x-scale
Frequency_MG2_PRECTFREQPDF=np.mean(P3_array,axis=1)*bounds_log10diff/bounds_lindiff
Frequency_P3_PRECTFREQPDF=np.mean(P3_array,axis=1)*bounds_log10diff/bounds_lindiff
Frequency_TRMM_PRECTFREQPDF=np.mean(TRMM_array,axis=1)*bounds_log10diff/bounds_lindiff

fig=plt.figure(figsize=(7,6))
ax=fig.add_subplot(1,1,1)
plt.plot(bincenter,Frequency_MG2_PRECTFREQPDF*100.,color='tab:blue',linewidth=2,label='MG2')
plt.plot(bincenter,Frequency_P3_PRECTFREQPDF*100.,color='tab:orange',linewidth=2,label='P3')
plt.plot(bincenter,Frequency_TRMM_PRECTFREQPDF*100.,color='tab:brown',linewidth=2,label='TRMM')
plt.ylabel('Frequency (%)',fontsize=16)
plt.xlabel('P (mm/day)',fontsize=16)
plt.yticks(size=16)
plt.xticks(size=16)
plt.yscale('log') 
plt.xlim((1,600))
plt.title('Frequency PDF format 30N-30S',fontsize=16)
plt.legend(loc='upper center',fontsize=16)
plt.savefig('Figurename')   #user modify

# =========================

# Do the same as above, but for CONUS box

# ==========================


# Retrieve Tropical mean PDF of FREQ and AMNT PDFs                                                                                
bincenter,MG2_2011_conus_PRECTFREQPDF,MG2_2011_conus_PRECTAMNTPDF=WC_diag_amwg.CONUS_PDF_values(f_in_MG2_2011,'PRECT')
bincenter,MG2_2012_conus_PRECTFREQPDF,MG2_2012_conus_PRECTAMNTPDF=WC_diag_amwg.CONUS_PDF_values(f_in_MG2_2012,'PRECT')
bincenter,MG2_2013_conus_PRECTFREQPDF,MG2_2013_conus_PRECTAMNTPDF=WC_diag_amwg.CONUS_PDF_values(f_in_MG2_2013,'PRECT')


bincenter,TRMM_2000_conus_PRECTFREQPDF,TRMM_2000_conus_PRECTAMNTPDF  =WC_diag_amwg.CONUS_PDF_values(f_in_TRMM_PDF_2000,'PRECT')
bincenter,TRMM_2001_conus_PRECTFREQPDF,TRMM_2001_conus_PRECTAMNTPDF  =WC_diag_amwg.CONUS_PDF_values(f_in_TRMM_PDF_2001,'PRECT')
bincenter,TRMM_2002_conus_PRECTFREQPDF,TRMM_2002_conus_PRECTAMNTPDF  =WC_diag_amwg.CONUS_PDF_values(f_in_TRMM_PDF_2002,'PRECT')


MG2_array=np.zeros((129,3))
MG2_array[:,0]=MG2_2011_conus_PRECTFREQPDF
MG2_array[:,1]=MG2_2012_conus_PRECTFREQPDF
MG2_array[:,2]=MG2_2013_conus_PRECTFREQPDF


TRMM_array=np.zeros((129,3))
TRMM_array[:,0]=TRMM_2000_conus_PRECTFREQPDF
TRMM_array[:,1]=TRMM_2001_conus_PRECTFREQPDF
TRMM_array[:,2]=TRMM_2002_conus_PRECTFREQPDF

Frequency_MG2_CONUS_PRECTFREQPDF=np.mean(MG2_array,axis=1)*bounds_log10diff/bounds_lindiff
Frequency_TRMM_CONUS_PRECTFREQPDF=np.mean(TRMM_array,axis=1)*bounds_log10diff/bounds_lindiff


fig=plt.figure(figsize=(7,6))
ax=fig.add_subplot(1,1,1)
plt.plot(bincenter,Frequency_MG2_CONUS_PRECTFREQPDF*100.,color='tab:blue',linewidth=2,label='MG2')    #user modify
plt.plot(bincenter,Frequency_TRMM_CONUS_PRECTFREQPDF*100.,color='tab:brown',linewidth=2,label='TRMM')
plt.ylabel('Frequency (%)',fontsize=16)
plt.xlabel('P (mm/day)',fontsize=16)
plt.yticks(size=16)
plt.xticks(size=16)
plt.yscale('log') 
plt.xlim((1,210))
plt.title('Frequency PDF format Continental US (35-49N, 125-75W) (3-hourly)',fontsize=16)
plt.legend(loc='upper center',fontsize=16)
plt.savefig('Figurename')   #user modify
