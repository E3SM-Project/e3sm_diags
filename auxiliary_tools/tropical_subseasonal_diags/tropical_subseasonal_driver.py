from __future__ import annotations

import glob
import json
import os
from typing import TYPE_CHECKING  # , Optional

import numpy as np
import xarray as xr

import e3sm_diags
from e3sm_diags.driver import utils
from e3sm_diags.logger import custom_logger
from e3sm_diags.plot.cartopy.mp_partition_plot import plot
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from zwf import zwf_functions as wf

logger = custom_logger(__name__)

# Script to compute and plot spectral powers of a subseasonal tropical field in 
#   zonal wavenumber-frequency space.  Both the plot files and files containing the
#   associated numerical data shown in the plots are created.

# Authors: Jim Benedict and Brian Medeiros
# Modified by Jill Zhang to integrate into E3SM Diags. 

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]
    """Return index of [array] closest in value to [value]
    Example:
    array = [ 0.21069679  0.61290182  0.63425412  0.84635244  0.91599191  0.00213826
              0.17104965  0.56874386  0.57319379  0.28719469]
    print(find_nearest(array, value=0.5))
    # 0.568743859261

    """

def wf_analysis(x, **kwargs):
    """Return zonal wavenumber-frequency power spectra of x.  The returned spectra are:
    spec_sym:    Raw (non-normalized) power spectrum of the component of x that is symmetric about the equator.
    spec_asym:   Raw (non-normalized) power spectrum of the component of x that is antisymmetric about the equator.
    nspec_sym:   Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is symmetric about the equator.
    nspec_asym:  Normalized (by a smoothed red-noise background spectrum) power spectrum of the component of x that is antisymmetric about the equator.
    
    The NCL version of 'wkSpaceTime' smooths the symmetric and antisymmetric components
    along the frequency dimension using a 1-2-1 filter once.  
    
    """
    # Get the "raw" spectral power
    # OPTIONAL kwargs: 
    # segsize, noverlap, spd, latitude_bounds (tuple: (south, north)), dosymmetries, rmvLowFrq

    z2 = wf.spacetime_power(x, **kwargs)
    z2avg = z2.mean(dim='component')
    z2.loc[{'frequency':0}] = np.nan # get rid of spurious power at \nu = 0 (mean)
   
    # Following NCL's wkSpaceTime, apply one pass of a 1-2-1 filter along the frequency
    #   domain to the raw (non-normalized) spectra/um.
    #   Do not use 0 frequency when smoothing here.
    #   Use weights that sum to 1 to ensure that smoothing is conservative.
    z2s = wf.smoothFrq121(z2,1)

    # The background is supposed to be derived from both symmetric & antisymmetric
    # Inputs to the background spectrum calculation should be z2avg
    background = wf.smoothBackground_wavefreq(z2avg)
    # separate components
    spec_sym = z2s[0,...]
    spec_asy = z2s[1,...]
    # normalize:  Following NCL's wkSpaceTime, use lightly smoothed version of spectra/um
    #             as numerator
    nspec_sym = spec_sym / background
    nspec_asy = spec_asy / background

    #spec_sym.rename('spec_raw_sym')
    #spec_asy.rename('spec_raw_asy')
    #nspec_sym.rename('spec_norm_sym')
    #nspec_asy.rename('spec_norm_asy')
    #background.rename('spec_background')
    #return spec_sym, spec_asy, nspec_sym, nspec_asy, background
    return spec_sym.rename('spec_raw_sym'), spec_asy.rename('spec_raw_asy'), nspec_sym.rename('spec_norm_sym'), nspec_asy.rename('spec_norm_asy'), background.rename('spec_background')

def plot_spectrum(s, ofil=None, clevs=None, cmapSpec='viridis', component=None, spec_type = None, 
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of non-normalized (raw) symmetric power spectrum with shallow water curves."""

    PlotDesc = {}
    PlotDesc['spec_raw_sym'] = {"long_name_desc": f"{varName}: log-base10 of lightly smoothed spectral power of component symmetric about equator", "ref_fig_num": "Figure 1"} # Figure number from  Wheeler and Kiladis (1999)
    PlotDesc['spec_raw_asy'] = {"long_name_desc": f"{varName}: log-base10 of lightly smoothed spectral power of component antisymmetric about equator", "ref_fig_num": "Figure 1"} 
    PlotDesc['spec_norm_sym'] = {"long_name_desc": f"{varName}: lightly smoothed spectral power of component symmetric about equator, normalized by heavily smoothed background spectrum", "ref_fig_num": "Figure 3"} 
    PlotDesc['spec_norm_asy'] = {"long_name_desc": f"{varName}: lightly smoothed spectral power of component antisymmetric about equator, normalized by heavily smoothed background spectrum", "ref_fig_num": "Figure 3"} 
    PlotDesc['spec_background'] = {"long_name_desc": f"{varName}: heavily smoothed version of the mean of spectral powers associated with the components symmetric and antisymmetric about equator", "ref_fig_num": "Figure 2"} 
    print(s)
    print('xxxxx', s.name, PlotDesc[s.name])
    
    text_offset = 0.005
    fb  = [0, .8]    # frequency bounds for plot
    wnb = [-15, 15]  # zonal wavenumber bounds for plot
    if(max(s['frequency'].values) == 0.5):
      fb = [0, .5]
    if(do_zoom):
      fb  = [0, .18]
      wnb = [-7, 7]
    # get data for dispersion curves:
    swfreq,swwn = wf.genDispersionCurves(Ahe=equivDepths)
    # swfreq.shape # -->(6, 3, 50)
    if spec_type == "normalized": # with dispersion curves to plot
        # For n=1 ER waves, allow dispersion curves to touch 0 -- this is for plot aesthetics only
        for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for ER waves
          indMinPosFrqER = np.where(swwn[3,i,:] >= 0., swwn[3,i,:], 1e20).argmin()  # index of swwn for least positive wn
          swwn[3,i,indMinPosFrqER],swfreq[3,i,indMinPosFrqER] = 0.,0.    # this sets ER's frequencies to 0. at wavenumber 0.

    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)

    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)

    # Final data refinement:  transpose and trim, set 0 freq to NaN, take log10 for raw, refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(*wnb))
    z.loc[{'frequency':0}] = np.nan

    if 'spec_raw' in s.name: 

        east_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(1,3)).sum()
        west_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(-3,-1)).sum()
        ew_ratio   = east_power / west_power
        print("\neast_power: %12.5f" % east_power)
        print("west_power: %12.5f" % west_power)
        print("ew_ratio: %12.5f\n" % ew_ratio)
        
        z = np.log10(z)

    z.attrs["long_name"] = PlotDesc[s.name]["long_name_desc"]
    z.attrs["method"] = f"Follows {PlotDesc[s.name]['ref_fig_num']} methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"

    if 'spec_raw' in s.name:

        z.attrs["ew_ratio_method"] = "Sum of raw (not log10) symmetric spectral power for ZWNs +/- 1-3, periods 24-96 days"
        z.attrs["east_power"] = east_power.values
        z.attrs["west_power"] = west_power.values
        z.attrs["ew_ratio"]   = ew_ratio.values

    # Save plotted data z to file as xArray data array
    dataDesc = f"spectral_power_{srcID}_{vari}_{spec_type}_{component}_200101_201412"
    z.to_netcdf(outDataDir + "/"+ dataDesc + ".nc")

    fig, ax = plt.subplots()
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
    #img = ax.contourf(kmesh0, vmesh0, z, levels=np.linspace(0.2, 3.0, 16), cmap='Spectral_r',  extend='both')
    img = ax.contourf(kmesh0, vmesh0, z, levels=clevs, cmap=cmapSpecUse,  norm=normSpecUse, extend='both')
    img2 = ax.contour(kmesh0, vmesh0, z, levels=clevs, linewidths=1., linestyles='solid', colors='gray', alpha=0.7)
    ax.axvline(0, linestyle='dashed', color=perfrq_col, linewidth=perfrq_thk, alpha=disp_alpha)
    if( (1./30.) < fb[1] ):
      ax.axhline((1./30.), linestyle='dashed', color=perfrq_col, alpha=perfrq_alpha)
      ax.text(wnb[0]+1,(1./30.)+text_offset,'30 days',color=perfrq_col, alpha=perfrq_alpha)
    if( (1./6.) < fb[1] ):
      ax.axhline((1./6.), linestyle='dashed', color=perfrq_col, alpha=perfrq_alpha)
      ax.text(wnb[0]+1,(1./6.)+text_offset,'6 days',color=perfrq_col, alpha=perfrq_alpha)
    if( (1./3.) < fb[1] ):
      ax.axhline((1./3.), linestyle='dashed', color=perfrq_col, alpha=perfrq_alpha)
      ax.text(wnb[0]+1,(1./3.)+text_offset,'3 days',color=perfrq_col, alpha=perfrq_alpha)
    for ii in range(3,6):
        ax.plot(swk[ii, 0,:], swf[ii,0,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 1,:], swf[ii,1,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 2,:], swf[ii,2,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
    ax.set_xlim(wnb)
    ax.set_ylim(fb)
    #ax.set_title(varName + ": Log ${\sum_{15^{\circ}S}^{15^{\circ}N} Power_{SYM}}$")   # Version w/ LaTeX
    if spec_type == "raw":
        ax.set_title(f"{varName}: Log{{Sum(Power) from 15°S-15°N}}\n")       # Version w/o LaTeX
    elif spec_type == "normalized":
        ax.set_title(f"{varName}: {{Sum(Power) from 15°S-15°N}}/Background\n")          # Version w/o LaTeX
    else:
        ax.set_title(f"{varName}: Log{{Smoothed Background Power}}\n")

    ax.set_title(sourceID, loc='left')
    ax.set_title(f"{component}", loc='right')

    if spec_type == "normalized":
        # Shallow water dispersion curve line labels:  See https://matplotlib.org/stable/tutorials/text/text_intro.html
        # n=1 ER dispersion curve labels
        text_opt = {'fontsize': 9,'verticalalignment': 'center','horizontalalignment': 'center','clip_on': True,'bbox': {'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0}}
        iwave, ih = 3, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -11.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 3, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -9.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 3, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        ax.text(-7.,0.10,'n=1 ER',text_opt)

        # Kelvin dispersion curve labels
        iwave, ih = 4, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 4, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 10.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 4, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 14.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        ax.text(6.,0.13,'Kelvin',text_opt)

        # IG dispersion curve labels
        iwave, ih = 5, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 5, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        iwave, ih = 5, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',text_opt)
        ax.text(-10.,0.48,'n=1 WIG',text_opt)
        ax.text(5.,0.48,'n=1 EIG',text_opt)

        # MJO label
        ax.text(6.,0.0333,'MJO',text_opt)


    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    # Save fig
    fig.savefig(outDataDir + "/"+ dataDesc + "_plot.png", bbox_inches='tight', dpi=300)

    print("Plot file created: %s\n" % outDataDir + "/"+ dataDesc + "_plot.png")



#
# Options ... right now these only go into wk.spacetime_power()
#
do_zooming = False    # Set to True to also make plots to zoom into MJO spectral region,
                     #   in addition to the default (larger) spectral region
latBound = (-15,15)  # latitude bounds for analysis
spd      = 1    # SAMPLES PER DAY
nDayWin  = 96   # Wheeler-Kiladis [WK] temporal window length (days)
nDaySkip = -60  # time (days) between temporal windows [segments]
                # negative means there will be overlapping temporal segments
twoMonthOverlap = -1*nDaySkip

vari = "PRECT"
srcID = "model"
outDataDir = "/global/cfs/cdirs/e3sm/www/chengzhu/tests/tropical_diags"

opt = {'segsize': nDayWin,
       'noverlap': twoMonthOverlap,
       'spd': spd,
       'latitude_bounds': latBound,
       'dosymmetries': True,
       'rmvLowFrq':True}

# Generic settings for spectra plots -- do not need to ba adjusted by user, but can be
#   customized.  The user can add additional key-value pairs as necessary to expand
#   plotting customization.
contour_levs_raw_spec  = (-1.4,-1.2,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2)
#contour_levs_norm_spec = (0.9,1,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7,3,3.5,4)
cmapSpecR = ["white",
             "paleturquoise","lightblue","skyblue",
             "lightgreen","limegreen","green","darkgreen",
             "yellow","orange","orangered",
             "red","maroon","magenta","orchid","pink",
             "lavenderblush"]
contour_levs_norm_spec = (0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7,3.0)

cmapSpecN = ["white",
             "gainsboro","lightgray","silver",
             "paleturquoise","skyblue",
             "lightgreen","mediumseagreen","seagreen",
             "yellow","orange",
             "red","maroon","pink"]

optPlot = {'varName': vari,
           'sourceID': srcID,
           'do_zoom' : False,
           'disp_col': 'black',
           'disp_thk': 1.5,
           'disp_alpha': 0.60,
           'perfrq_col': 'dimgray',
           'perfrq_thk': 1.0,
           'perfrq_alpha': 0.80,
           'equivDepths' : [50, 25, 12]}



datapath = '/global/cfs/cdirs/e3sm/forsyth/E3SMv2/v2.LR.historical_0201/post/atm/180x360_aave/ts/daily/5yr'
data = xr.open_mfdataset(glob.glob(f"{datapath}/PRECT_201001_201412.nc")).sel(
        lat=slice(-15, 15))['PRECT']
# TODO: subset time

# Unit conversion
if vari == "PRECT":
  if data.attrs["units"] == "m/s" or data.attrs["units"] == "m s{-1}":
    print("\nBEFORE unit conversion: Max/min of data: " + str(data.values.max()) + "   " + str(data.values.min()))
    data.values = data.values * 1000. * 86400.    # convert m/s to mm/d, do not alter metadata (yet)
    data.attrs["units"] = "mm/d"                  # adjust metadata to reflect change in units
    print("\nAFTER unit conversion: Max/min of data: " + str(data.values.max()) + "   " + str(data.values.min()))
    print(data)
if vari == "precipAvg":
  if data.attrs["units"] == "mm/hr":
    print("\nBEFORE unit conversion: Max/min of data: " + str(data.values.max()) + "   " + str(data.values.min()))
    data.values = data.values * 24.               # convert mm/hr to mm/d, do not alter metadata (yet)
    data.attrs["units"] = "mm/d"                  # adjust metadata to reflect change in units
    print("\nAFTER unit conversion: Max/min of data: " + str(data.values.max()) + "   " + str(data.values.min()))
    print(data)


# Wavenumber Frequency Analysis
spec_raw_sym, spec_raw_asym, spec_norm_sym, spec_norm_asym, spec_background = wf_analysis(data, **opt)
print("\nspec_raw_sym metadata:")
print(spec_raw_sym.dims)
print(spec_raw_sym.coords)
print(spec_raw_sym.attrs)
print(spec_raw_sym.max())


outPlotName = outDataDir
plot_spectrum(spec_raw_sym, outPlotName, contour_levs_raw_spec, cmapSpecR,component="symmetric", spec_type = "raw", **optPlot)
plot_spectrum(spec_raw_asym, outPlotName,  contour_levs_raw_spec, cmapSpecR,component="antisymmetric", spec_type = "raw", **optPlot)
plot_spectrum(spec_background, outPlotName,  contour_levs_raw_spec, cmapSpecR,component="background", spec_type = "background", **optPlot)
plot_spectrum(spec_norm_sym, outPlotName, contour_levs_norm_spec, cmapSpecN,component="symmetric", spec_type = "normalized", **optPlot)
plot_spectrum(spec_norm_asym, outPlotName, contour_levs_norm_spec, cmapSpecN,component="antisymmetric", spec_type = "normalized", **optPlot)



