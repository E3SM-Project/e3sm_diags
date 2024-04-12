# Script to compute and plot spectral powers of a subseasonal tropical field in 
#   zonal wavenumber-frequency space.  Both the plot files and files containing the
#   associated numerical data shown in the plots are created.
#
# User-defined inputs can be entered near the beginning of the main script (the line
#   if __name__ == "__main__":)
#
# To invoke on NERSC-Cori or LCRC:
#   First, activate E3SM unified environment.  E.g, for NERSC-Cori:
#   > source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_cori-knl.sh
#   Then, simply run:
#   > python zwf_plot.py

import sys
import os
import math
import glob
import numpy as np
import xarray as xr
# our local module:
import zwf_functions as wf
#import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm


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
    return spec_sym, spec_asy, nspec_sym, nspec_asy, background


def plot_raw_symmetric_spectrum(s, ofil=None, dataDesc=None, clevs=None, cmapSpec='viridis', 
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of non-normalized (raw) symmetric power spectrum with shallow water curves."""

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
    
    # Final data refinement:  transpose and trim, set 0 freq to NaN, take log10, refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(*wnb))
    z.loc[{'frequency':0}] = np.nan
    east_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(1,3)).sum()
    west_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(-3,-1)).sum()
    ew_ratio   = east_power / west_power
    print("\neast_power: %12.5f" % east_power)
    print("west_power: %12.5f" % west_power)
    print("ew_ratio: %12.5f\n" % ew_ratio)
    
    z = np.log10(z)
    z.attrs["long_name"] = varName + ": log-base10 of lightly smoothed spectral power of component symmetric about equator"
    z.attrs["method"] = "Follows Figure 1 methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"
    z.attrs["ew_ratio_method"] = "Sum of raw (not log10) symmetric spectral power for ZWNs +/- 1-3, periods 24-96 days"
    z.attrs["east_power"] = east_power.values
    z.attrs["west_power"] = west_power.values
    z.attrs["ew_ratio"]   = ew_ratio.values
    
    
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
    ax.set_title(f"{varName}: Log{{Sum(Power) from 15°S-15°N}}\n")       # Version w/o LaTeX
    ax.set_title(sourceID, loc='left')
    ax.set_title("Symmetric", loc='right')

    if(not do_zoom):         # For now, only add equivalent depth and shallow water curve labels -NOT- to zoomed-in plots
        # Shallow water dispersion curve line labels:  See https://matplotlib.org/stable/tutorials/text/text_intro.html
        # n=1 ER dispersion curve labels
        iwave, ih = 3, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -11.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 3, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -9.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 3, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-7.,0.10,'n=1 ER',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
     
        # Kelvin dispersion curve labels
        iwave, ih = 4, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 4, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 10.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 4, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 14.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(6.,0.13,'Kelvin',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # IG dispersion curve labels
        iwave, ih = 5, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 5, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 5, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-10.,0.48,'n=1 WIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(5.,0.48,'n=1 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # MJO label
        ax.text(6.,0.0333,'MJO',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})

    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    if ofil is not None:
        fig.savefig(ofil, bbox_inches='tight', dpi=300)
        print("Plot file created: %s\n" % ofil)
    
    # Save plotted data z to file as xArray data array
    if(not do_zoom):
      z.to_netcdf(outDataDir + "/zwfData_raw_sym_" + dataDesc + ".nc")
    

def plot_normalized_symmetric_spectrum(s, ofil=None, dataDesc=None, clevs=None, cmapSpec='Spectral_r', 
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of normalized symmetric power spectrum with shallow water curves."""
    
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
    # For n=1 ER waves, allow dispersion curves to touch 0 -- this is for plot aesthetics only
    for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for ER waves
      indMinPosFrqER = np.where(swwn[3,i,:] >= 0., swwn[3,i,:], 1e20).argmin()  # index of swwn for least positive wn
      swwn[3,i,indMinPosFrqER],swfreq[3,i,indMinPosFrqER] = 0.,0.    # this sets ER's frequencies to 0. at wavenumber 0.
    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)
    
#    for i in range(6):
#      print(f'\nwaveType={i}')
#      print(f"{'waveNum':12}  {'frq (ih=0)':12}  {'frq (ih=1)':12}  {'frq (ih=2)':12}")
#      for j in range(50):
#        print(f'{swk[i,0,j]:12.4f}  {swf[i,0,j]:12.4f}  {swf[i,1,j]:12.4f}  {swf[i,2,j]:12.4f}')
#    sys.exit()
        
    
    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)
    
    # Final data refinement:  transpose and trim, set 0 freq to NaN (no log10 for
    #   normalized results), refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(-15,15))
    z.loc[{'frequency':0}] = np.nan
    z.attrs["long_name"] = varName + ": lightly smoothed spectral power of component symmetric about equator, normalized by heavily smoothed background spectrum"
    z.attrs["method"] = "Follows Figure 3 methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"
    
    fig, ax = plt.subplots()
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
    #img = ax.contourf(kmesh0, vmesh0, z, levels=np.linspace(0.2, 3.0, 16), cmap='Spectral_r',  extend='both')
    img = ax.contourf(kmesh0, vmesh0, z, levels=clevs, cmap=cmapSpecUse,  norm=normSpecUse,  extend='both')
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
    #ax.set_title(varName + ": $\sum_{15^{\circ}S}^{15^{\circ}N} Power_{SYM}$ / Background")       # Version w/ LaTeX
    ax.set_title(f"{varName}: {{Sum(Power) from 15°S-15°N}}/Background\n")          # Version w/o LaTeX
    ax.set_title(sourceID, loc='left')
    ax.set_title("Symmetric", loc='right')
    
    if(not do_zoom):         # For now, only add equivalent depth and shallow water curve labels -NOT- to zoomed-in plots
        # Shallow water dispersion curve line labels:  See https://matplotlib.org/stable/tutorials/text/text_intro.html
        # n=1 ER dispersion curve labels
        iwave, ih = 3, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -11.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 3, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -9.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 3, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
       ax.text(-7.,0.10,'n=1 ER',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
     
        # Kelvin dispersion curve labels
        iwave, ih = 4, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 4, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 10.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 4, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 14.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(6.,0.13,'Kelvin',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # IG dispersion curve labels
        iwave, ih = 5, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 5, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 5, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 0.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-10.,0.48,'n=1 WIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(5.,0.48,'n=1 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # MJO label
        ax.text(6.,0.0333,'MJO',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    if ofil is not None:
        fig.savefig(ofil, bbox_inches='tight', dpi=300)
        print("Plot file created: %s\n" % ofil)
    
    # Save plotted data z to file as xArray data array
    if(not do_zoom):
      z.to_netcdf(outDataDir + "/zwfData_norm_sym_" + dataDesc + ".nc")
    

def plot_raw_asymmetric_spectrum(s, ofil=None, dataDesc=None, clevs=None, cmapSpec='viridis',
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of non-normalized (raw) antisymmetric power spectrum with shallow water curves."""
    
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
    # "Connect" MRG and n=0 IG dispersion curves across wavenumber 0 -- this is for plot aesthetics only
    for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for MRG and n=0 IG curves
      indPosWNclosest0 = np.where(swwn[0,i,:] >= 0., swwn[0,i,:], 1e20).argmin()  # index of swwn for positive wn closest to 0
      swfreq[0,i,indPosWNclosest0] = swfreq[1,i,indPosWNclosest0]    # this sets MRG's frequencies at least positive wn to n=0 IG's frequencies at least positive wn
    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)
    
    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)
    
    # Final data refinement:  transpose and trim, set 0 freq to NaN, take log10, refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(-15,15))
    z.loc[{'frequency':0}] = np.nan
    z = np.log10(z)
    z.attrs["long_name"] = varName + ": log-base10 of lightly smoothed spectral power of component antisymmetric about equator"
    z.attrs["method"] = "Follows Figure 1 methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"
    
    fig, ax = plt.subplots()
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
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
    for ii in range(0,3):
        ax.plot(swk[ii, 0,:], swf[ii,0,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 1,:], swf[ii,1,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 2,:], swf[ii,2,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
    ax.set_xlim(wnb)
    ax.set_ylim(fb)
    #ax.set_title(varName + ": Log10 $\sum_{15^{\circ}S}^{15^{\circ}N} Power_{ASYM}$")       # Version w/ LaTeX
    ax.set_title(f"{varName}: Log{{Sum(Power) from 15°S-15°N}}\n")       # Version w/o LaTeX
    ax.set_title(sourceID, loc='left')
    ax.set_title("Antisymmetric", loc='right')

    if(not do_zoom):
        # Shallow water dispersion curve line labels:  See https://matplotlib.org/stable/tutorials/text/text_intro.html
        # MRG dispersion curve labels -- SKIP LABELING EQUIVALENT DEPTHS FOR MRG WAVES AND ONLY LABEL FOR N=2 EIG WAVES, WHICH ARE POSTIVE-WAVENUMBER EXTENSIONS OF THE MRG CURVES
#        iwave, ih = 0, 0
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 4.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
#        iwave, ih = 0, 1
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 6.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
#        iwave, ih = 0, 2
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-6.,0.18,'MRG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # n=0 EIG dispersion curve labels
        iwave, ih = 1, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 5.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 1, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 1, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(9.,0.48,'n=0 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # n=2 IG dispersion curve labels
        iwave, ih = 2, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 2, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 2, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-10.,0.65,'n=2 WIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(8.,0.65,'n=2 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # MJO label
        ax.text(3.,0.0333,'MJO',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})

    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    if ofil is not None:
        fig.savefig(ofil, bbox_inches='tight', dpi=300)
        print("Plot file created: %s\n" % ofil)
    
    # Save plotted data z to file as xArray data array
    if(not do_zoom):
      z.to_netcdf(outDataDir + "/zwfData_raw_asym_" + dataDesc + ".nc")
    

def plot_normalized_asymmetric_spectrum(s, ofil=None, dataDesc=None, clevs=None, cmapSpec='Spectral_r',
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of normalized antisymmetric power spectrum with shallow water curves."""
    
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
    # "Connect" MRG and n=0 IG dispersion curves across wavenumber 0 -- this is for plot aesthetics only
    for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for MRG and n=0 IG curves
      indPosWNclosest0 = np.where(swwn[0,i,:] >= 0., swwn[0,i,:], 1e20).argmin()  # index of swwn for positive wn closest to 0
      swfreq[0,i,indPosWNclosest0] = swfreq[1,i,indPosWNclosest0]    # this sets MRG's frequencies at least positive wn to n=0 IG's frequencies at least positive wn
    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)
    
    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)
    
    # Final data refinement:  transpose and trim, set 0 freq to NaN (no log10 for
    #   normalized results), refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(-15,15))
    z.loc[{'frequency':0}] = np.nan
    z.attrs["long_name"] = varName + ": lightly smoothed spectral power of component antisymmetric about equator, normalized by heavily smoothed background spectrum"
    z.attrs["method"] = "Follows Figure 3 methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"
        
    fig, ax = plt.subplots()
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
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
    for ii in range(0,3):
        ax.plot(swk[ii, 0,:], swf[ii,0,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 1,:], swf[ii,1,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
        ax.plot(swk[ii, 2,:], swf[ii,2,:], color=disp_col, linewidth=disp_thk, alpha=perfrq_alpha)
    ax.set_xlim(wnb)
    ax.set_ylim(fb)    
    #ax.set_title(varName + ": $\sum_{15^{\circ}S}^{15^{\circ}N} Power_{SYM}$ / Background")       # Version w/ LaTeX
    ax.set_title(f"{varName}: {{Sum(Power) from 15°S-15°N}}/Background\n")          # Version w/o LaTeX
    ax.set_title(sourceID, loc='left')
    ax.set_title("Antisymmetric", loc='right')
    
    if(not do_zoom):
        # Shallow water dispersion curve line labels:  See https://matplotlib.org/stable/tutorials/text/text_intro.html
        # MRG dispersion curve labels -- SKIP LABELING EQUIVALENT DEPTHS FOR MRG WAVES AND ONLY LABEL FOR N=2 EIG WAVES, WHICH ARE POSTIVE-WAVENUMBER EXTENSIONS OF THE MRG CURVES
#        iwave, ih = 0, 0
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 4.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
#        iwave, ih = 0, 1
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 6.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
#        iwave, ih = 0, 2
#        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
#        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-6.,0.18,'MRG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # n=0 EIG dispersion curve labels
        iwave, ih = 1, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 5.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 1, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 1, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], 8.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(9.,0.48,'n=0 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # n=2 IG dispersion curve labels
        iwave, ih = 2, 0
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 2, 1
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        iwave, ih = 2, 2
        idxClose,valClose = find_nearest(swk[iwave,ih,:], -2.)    # Locate index of wavenumber closest to input value [and the actual (float) wavenumber value]
        ax.text(valClose,swf[iwave,ih,idxClose],f'{equivDepths[ih]}',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(-10.,0.65,'n=2 WIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
        ax.text(8.,0.65,'n=2 EIG',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
        # MJO label
        ax.text(3.,0.0333,'MJO',fontsize=9,verticalalignment='center',horizontalalignment='center',clip_on=True,bbox={'facecolor': 'white', 'edgecolor':'none', 'alpha': 0.7, 'pad': 0.0})
    
    
    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    if ofil is not None:
        fig.savefig(ofil, bbox_inches='tight', dpi=300)
        print("Plot file created: %s\n" % ofil)
    
    # Save plotted data z to file as xArray data array
    if(not do_zoom):
      z.to_netcdf(outDataDir + "/zwfData_norm_asym_" + dataDesc + ".nc")
    

def plot_background_spectrum(s, ofil=None, dataDesc=None, clevs=None, cmapSpec='viridis', 
        varName=None, sourceID=None, do_zoom=False,
        disp_col='black', disp_thk=1.5, disp_alpha=0.60,
        perfrq_col='dimgray', perfrq_thk=1.0, perfrq_alpha=0.80, equivDepths=[50, 25, 12]):
    """Basic plot of background power spectrum with shallow water curves."""
    
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
    # For n=1 ER waves, allow dispersion curves to touch 0 -- this is for plot aesthetics only
    for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for ER waves
      indMinPosFrqER = np.where(swwn[3,i,:] >= 0., swwn[3,i,:], 1e20).argmin()  # index of swwn for least positive wn
      swwn[3,i,indMinPosFrqER],swfreq[3,i,indMinPosFrqER] = 0.,0.    # this sets ER's frequencies to 0. at wavenumber 0.
    # "Connect" MRG and n=0 IG dispersion curves across wavenumber 0 -- this is for plot aesthetics only
    for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for MRG and n=0 IG curves
      indPosWNclosest0 = np.where(swwn[0,i,:] >= 0., swwn[0,i,:], 1e20).argmin()  # index of swwn for positive wn closest to 0
      swfreq[0,i,indPosWNclosest0] = swfreq[1,i,indPosWNclosest0]    # this sets MRG's frequencies at least positive wn to n=0 IG's frequencies at least positive wn
    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)
    
    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)
    
    # Final data refinement:  transpose and trim, set 0 freq to NaN, take log10, refine metadata
    z  = s.transpose().sel(frequency=slice(*fb), wavenumber=slice(-15,15))
    z.loc[{'frequency':0}] = np.nan
    z = np.log10(z)
    z.attrs["long_name"] = varName + ": heavily smoothed version of the mean of spectral powers associated with the components symmetric and antisymmetric about equator"
    z.attrs["method"] = "Follows Figure 2 methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"
    
    fig, ax = plt.subplots()
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
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
    ax.set_xlim(wnb)
    ax.set_ylim(fb)
    ax.set_title(f"{varName}: Log{{Smoothed Background Power}}\n")
    ax.set_title(sourceID, loc='left')
#    ax.set_title("", loc='right')
    plt.ylabel("Frequency (CPD)")
    plt.xlabel("Zonal wavenumber")
    plt.gcf().text(0.12, 0.03, "Westward", fontsize=11)
    plt.gcf().text(0.64, 0.03, "Eastward", fontsize=11)
    fig.colorbar(img)
    if ofil is not None:
        fig.savefig(ofil, bbox_inches='tight', dpi=300)
        print("Plot file created: %s\n" % ofil)

    # Save plotted data z to file as xArray data array
    if(not do_zoom):
      z.to_netcdf(outDataDir + "/zwfData_background_" + dataDesc + ".nc")

#
# LOAD DATA, x = DataArray(time, lat, lon), e.g., daily mean precipitation
#
def get_data(filenames, variablename):
    if(len(filenames) == 1):
      print("\nOpening " + variablename + " from single file: " + filenames[0])
      try: 
          ds = xr.open_dataset(filenames[0])
      except ValueError:
          ds = xr.open_dataset(filenames[0], decode_times=False)
    else:
      print("\nOpening " + variablename + " across multiple files: ")
      for j in range(len(filenames)):
         print("   " + str(j) + "  " + filenames[j])
      try: 
#          ds = xr.open_mfdataset(filenames, combine="by_coords")
          ds = xr.open_mfdataset(filenames)
      except:
          print("Error reading multi-file data set")
    x = ds[variablename].compute()
#    return ds[variablename]
    return x


if __name__ == "__main__":
#    swfreq,swwn = wf.genDispersionCurves()
#    for i in range(6):
#      print(f'\nwaveType={i}')
#      print(f"{'waveNum':12}  {'frq (ih=0)':12}  {'frq (ih=1)':12}  {'frq (ih=2)':12}")
#      for j in range(50):
#        print(f'{swwn[i,0,j]:12.4f}  {swfreq[i,0,j]:12.4f}  {swfreq[i,1,j]:12.4f}  {swfreq[i,2,j]:12.4f}')
#    sys.exit()
    #
    # ==================   User parameters (below)   ============================
    #
    # input file -- could make this a CLI argument
    #   * fili could specify a single (pre-concatenated) file, or a directory containing
    #     a collection of files.  If it's set to a directory, do -not- include trailing "/".
    #     If it's set to a directory, yrBeg and yrEnd will defined the
    #     (inclusive) year span from which to read multiple files and concatenate them
    #     along the time dimension.
    #   * vari:  The name of the variable to read in, from fili
    #   * yrBeg, yrEnd:  The beginning and ending years of the data
    #
    vari  = "precipAvg"  # "PRECT"
#    fili = "/lcrc/group/e3sm/ac.golaz/E3SMv2/v2.LR.piControl/archive/atm/hist"

#    srcID = "TRMM"
#    fili = "/global/cfs/cdirs/e3sm/benedict/obs/trmm_1dd.h2.20010101_to_20101231.PRECT.nc"
##    fili = "/lcrc/group/acme/ac.benedict/test_data_in/trmm_1dd.h2.20010101_to_20101231.PRECT.nc"
#    yrBeg = 2001
#    yrEnd = 2010
#    histStr = ""   # Input file descriptors, where file name is [histStr].yyyy-mm-dd*.nc
    
    srcID = "IMERG_2.5deg"
    fili = "/global/cfs/cdirs/e3sm/benedict/mjo_isv/obs_reanalysis/IMERG/daily_73x144"
    yrBeg = 2002
    yrEnd = 2006
    histStr = "3B-DAY.MS.MRG.3IMERG.V06"   # Input file descriptors, where file name is [histStr].yyyy-mm-dd*.nc

#    srcID = "E3SMv1_0.25deg"
#    fili = "/global/cscratch1/sd/ruplanl9/e3sm_scratch/highRes_1980_2015/PRECT/remap_to_0.25deg"
#    yrBeg = 2005
#    yrEnd = 2009
#    histStr = "202101027-maint-1.0-tro.A_WCYCL20TRS_CMIP6_HR.ne120_oRRS18v3_ICG.unc12.cam.h3"   # Input file descriptors, where file name is [histStr].yyyy-mm-dd*.nc

#    srcID = "DECK.v1b.LR.CMIP6"
#    fili = "/lcrc/group/acme/ac.benedict/test_data_in/20201125.DECKv1b_H1.cmip6_180x360.edison.daily.PRECT.nc"
#    yrBeg = 2005
#    yrEnd = 2014
#    histStr = ""   # Input file descriptors, where file name is [histStr].yyyy-mm-dd*.nc

#    srcID = "v2.LR.piControl"
#    fili = "/lcrc/group/acme/ac.benedict/test_data_in/remapped"
#    yrBeg = 491
#    yrEnd = 500
#    histStr = "v2.LR.piControl.eam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc

#    srcID = "v2.LR.historical_0101"
#    fili = "/lcrc/group/acme/ac.benedict/test_data_in/remapped"
#    yrBeg = 2005
#    yrEnd = 2014
#    histStr = "v2.LR.historical_0101.eam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc


### -------   cesm2.0.1   -------
#    srcID = "B1850_c201_CTL"
#    fili = "/global/cfs/cdirs/e3sm/benedict/cesm2/B1850_c201_CTL"
#    yrBeg = 1
#    yrEnd = 30
#    histStr = "B1850_c201_CTL.cam.h2"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc


### -------   v1   -------
#    srcID = "DECKv1b_H1"
#    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/20180215.DECKv1b_H1.ne30_oEC.edison/remapped"
#    yrBeg = 1985  # 2001
#    yrEnd = 2014  # 2010
#    histStr = "20180215.DECKv1b_H1.ne30_oEC.edison.cam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc
    
#    srcID = "DECKv1b_A1"
#    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/20180316.DECKv1b_A1.ne30_oEC.edison/remapped"
#    yrBeg = 1985  # 2001
#    yrEnd = 2014  # 2010
#    histStr = "20180316.DECKv1b_A1.ne30_oEC.edison.cam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc

#    srcID = "DECKv1b_piControl"
#    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/20180129.DECKv1b_piControl.ne30_oEC.edison/remapped"
#    yrBeg = 250
#    yrEnd = 279
#    histStr = "20180129.DECKv1b_piControl.ne30_oEC.edison.cam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc
    

### -------   v2   -------
#    srcID = "v2.LR.historical_0201"
#    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/v2.LR.historical_0201/remapped"
#    yrBeg = 1985  # 2001
#    yrEnd = 2014  # 2010
#    histStr = "v2.LR.historical_0201.eam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc
    
#    srcID = "v2.LR.amip_0201"
#    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/v2.LR.amip_0201/remapped"
#    yrBeg = 1985  # 2001
#    yrEnd = 2014  # 2010
#    histStr = "v2.LR.amip_0201.eam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc

#    srcID = "v2.LR.piControl"
#    fili = "/global/cfs/cdirs/e3sm/benedict/e3sm_v2_output/v2.LR.piControl/remapped"
##    fili = "/lcrc/group/acme/ac.benedict/e3sm_data/atm/v2.LR.piControl/remapped"
#    yrBeg = 460
#    yrEnd = 489
#    histStr = "PRECT_v2.LR.piControl.eam.h1"   # Input file descriptors, where file name is [vari]_[histStr].yyyy-mm-dd*.nc
    
    
    # outputs
    outDataDir = "/global/cfs/cdirs/e3sm/benedict/mjo_isv/data_E3SMdiags_TEST"
    outPlotDir = "/global/cfs/cdirs/e3sm/benedict/mjo_isv/plot_E3SMdiags_TEST"
#    outDataDir = "/home/ac.benedict/analysis/results/dataTrop"
#    outPlotDir = "/home/ac.benedict/analysis/results/plotTrop"

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
#    if(vari == "U850"):
#      contour_levs_raw_spec  = np.linspace(-2.6, 0.2, num=15)
#      contour_levs_norm_spec = (0.7,0.75,0.8,0.85,0.9,0.95,1,1.1,1.2,1.3,1.4,1.5,2,2.5,3)
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
    # disp_col      # Color for dispersion lines/labels
    # disp_thk      # Thickness for dispersion lines
    # disp_alpha    # Transparency for dispersion lines/labels (alpha=0 for opaque)
    # perfrq_col    # Color for period and frequency ref lines
    # perfrq_thk    # Thickness for period and frequency ref lines
    # perfrq_alpha  # Transparency for period and frequency ref lines (alpha=0 for opaque)
    # clevs         # Contour levels
    # ==================   User parameters (above)   ============================
    
    
    
    # Load data -- ORIGINAL SINGLE FILE
#    data = get_data(fili, vari)  # expected input is (time,lat,lon)
    
    # Load data:
    #   If fili is a single file, fileList will be a single-element list
    #   If fili is a directory, fileList will be a list containing multiple file names
    #     (contained within that directory).  The file names will span a year range
    #     bounded by yrBeg and yrEnd (inclusive).  ASSUMPTION:  The input data files
    #     have the standard YYYY-MM-DD date structure.
    if os.path.isfile(fili):
      fileList = [fili]
      data = get_data(fileList, vari)
    elif os.path.isdir(fili):
      yrList = [i for i in range(yrBeg,yrEnd+1)]
      print(yrList)
      yrListStr = []
      for y in yrList:
        yrListStr.append("%0.4i" % y)
      print(yrListStr)
      fileList = []
      for y in yrListStr:
#        fileList += sorted(glob.glob("{path}/*{hStr}.{year}*.nc".format(path=fili, hStr=histStr, year=y)))
        print("{path}/{hStr}.{year}*.nc".format(path=fili, hStr=histStr, year=y))
#        fileList += sorted(glob.glob("{path}/{hStr}.{year}-??-??*.nc".format(path=fili, hStr=histStr, year=y)))
        #fileList += sorted(glob.glob("{path}/{hStr}.{year}-??-??*.nc".format(path=fili, hStr=histStr, year=y)))
        fileList += sorted(glob.glob("{path}/{hStr}.{year}*.remap_73x144.nc".format(path=fili, hStr=histStr, year=y)))
#      print("\n\n\n JJB: fileList is:")
#      print(fileList)
      data = get_data(fileList, vari)
    else:
      sys.exit("\nFatal error: fili entry " + fili + " is neither a valid single file nor a valid directory. Exiting.")
    
    if not fileList:         # In python, empty lists are boolean 'false'
      sys.exit("\nFatal error: fileList is empty, no files matching the selected time range exist. Exiting.")
    
    print(data)
    
    # Assess/confirm time range of input data
    xt = data.coords['time']
    yrBegData = xt[0].time.dt.year
    yrEndData = xt[-1].time.dt.year
    if(xt[-1].time.dt.dayofyear < 183):     # if final time step is < halfway through year, subtract 1 from yrEndData
      yrEndData = yrEndData - 1
    print("\nyrBegData (may not represent full year!): %0.4d" % yrBegData)
    print("yrEndData (may not represent full year!): %0.4d" % yrEndData)
    
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
          
    opt = {'segsize': nDayWin, 
           'noverlap': twoMonthOverlap, 
           'spd': spd, 
           'latitude_bounds': latBound, 
           'dosymmetries': True, 
           'rmvLowFrq':True}    
    
    spec_raw_sym, spec_raw_asym, spec_norm_sym, spec_norm_asym, spec_background = wf_analysis(data, **opt)
    print("\nspec_raw_sym metadata:")
    print(spec_raw_sym.dims)
    print(spec_raw_sym.coords)
    print(spec_raw_sym.attrs)
    print(spec_raw_sym.max())
    
 
    #
    # Plots ... sort of matching NCL, but not worrying much about customizing.
    #
    # This string goes into the output file name encapsulating the data that is plotted:
    #   Note that 'yrBegData' and 'yrEndData' might not represent the full year
    plotFileDescriptor = "%s_%0.4d-%0.4d_%s" % (srcID, yrBegData, yrEndData, vari)
    
    # "Fig1": Log-10 of raw spectral power
    outPlotName = outPlotDir + "/" + "zwfPlot_raw_sym_" + plotFileDescriptor + ".png"
    plot_raw_symmetric_spectrum(spec_raw_sym, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
    outPlotName = outPlotDir + "/" + "zwfPlot_raw_asym_" + plotFileDescriptor + ".png"
    plot_raw_asymmetric_spectrum(spec_raw_asym, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
    # "Fig2": Log-10 of smoothed background spectrum
    outPlotName = outPlotDir + "/" + "zwfPlot_background_" + plotFileDescriptor + ".png"
    plot_background_spectrum(spec_background, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
    # "Fig3": Log-10 of normalized spectral power
    outPlotName = outPlotDir + "/" + "zwfPlot_norm_sym_" + plotFileDescriptor + ".png"
    plot_normalized_symmetric_spectrum(spec_norm_sym, outPlotName, plotFileDescriptor, contour_levs_norm_spec, cmapSpecN, **optPlot)
    
    outPlotName = outPlotDir + "/" + "zwfPlot_norm_asym_" + plotFileDescriptor + ".png"
    plot_normalized_asymmetric_spectrum(spec_norm_asym, outPlotName, plotFileDescriptor, contour_levs_norm_spec, cmapSpecN, **optPlot)

    
    # If also making plots to zoom into MJO spectral region
    if(do_zooming):
      
      optPlot['do_zoom'] = do_zooming    # Change do_zoom value in optPlot dictionary to True
      
      # "Fig1": Log-10 of raw spectral power
      outPlotName = outPlotDir + "/" + "zwfPlot_raw_sym_ZOOM_" + plotFileDescriptor + ".png"
      plot_raw_symmetric_spectrum(spec_raw_sym, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
      outPlotName = outPlotDir + "/" + "zwfPlot_raw_asym_ZOOM_" + plotFileDescriptor + ".png"
      plot_raw_asymmetric_spectrum(spec_raw_asym, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
      # "Fig2": Log-10 of smoothed background spectrum
      outPlotName = outPlotDir + "/" + "zwfPlot_background_ZOOM_" + plotFileDescriptor + ".png"
      plot_background_spectrum(spec_background, outPlotName, plotFileDescriptor, contour_levs_raw_spec, cmapSpecR, **optPlot)
    
      # "Fig3": Log-10 of normalized spectral power
      outPlotName = outPlotDir + "/" + "zwfPlot_norm_sym_ZOOM_" + plotFileDescriptor + ".png"
      plot_normalized_symmetric_spectrum(spec_norm_sym, outPlotName, plotFileDescriptor, contour_levs_norm_spec, cmapSpecN, **optPlot)
    
      outPlotName = outPlotDir + "/" + "zwfPlot_norm_asym_ZOOM_" + plotFileDescriptor + ".png"
      plot_normalized_asymmetric_spectrum(spec_norm_asym, outPlotName, plotFileDescriptor, contour_levs_norm_spec, cmapSpecN, **optPlot)
