from __future__ import annotations

import os
from typing import TYPE_CHECKING

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import xarray as xr
from matplotlib.colors import BoundaryNorm, ListedColormap

#from e3sm_diags.plot.utils import _add_colormap, _save_plot
from zwf import zwf_functions as wf

from e3sm_diags.logger import custom_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

if TYPE_CHECKING:
    from e3sm_diags.driver.lat_lon_driver import MetricsDict


matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = custom_logger(__name__)

# Plot title and side title configurations.
PLOT_TITLE = {"fontsize": 11.5}
PLOT_SIDE_TITLE = {"fontsize": 9.5}
        
# Position and sizes of subplot axes in page coordinates (0 to 1)
PANEL = [
    (0.17, 0.70, 0.50, 0.25),
    (0.17, 0.37, 0.50, 0.25),
    (0.17, 0.04, 0.50, 0.25),
]       
        
# Border padding relative to subplot axes for saving individual panels
# (left, bottom, right, top) in page coordinates
BORDER_PADDING = (-0.06, -0.03, 0.13, 0.03)

CONTOUR_LEVS_SPEC_RAW = (-1.4,-1.2,-1,-0.9,-0.8,-0.7,-0.6,-0.5,-0.4,-0.3,-0.2,-0.1,0,0.1,0.2)

CMAP_SPEC_RAW = ["white",
             "paleturquoise","lightblue","skyblue",
             "lightgreen","limegreen","green","darkgreen",
             "yellow","orange","orangered",
             "red","maroon","magenta","orchid","pink",
             "lavenderblush"]

CONTOUR_LEVS_SPEC_NORM = (0.9,1.0,1.1,1.2,1.3,1.4,1.5,1.7,1.9,2.1,2.4,2.7,3.0)
      
CMAP_SPEC_NORM = ["white",
             "gainsboro","lightgray","silver",
             "paleturquoise","skyblue",
             "lightgreen","mediumseagreen","seagreen",
             "yellow","orange",
             "red","maroon","pink"]

CONTOUR_LEVS_SPEC_RAW_DIFF = (-80.,-60.,-40.,-20.,-10.,-5.,5.,10.,20.,40.,60.,80.)
CONTOUR_LEVS_SPEC_NORM_DIFF = (-60.,-30.,-20.,-15.,-10.,-5.,5.,10.,15.,20.,30.,60.)


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

def create_colormap_clevs(cmapSpec, clevs):
    cmapSpecUse = ListedColormap(cmapSpec[1:-1])   # recall: range is NOT inclusive for final index
    cmapSpecUse.set_under(cmapSpec[0])
    cmapSpecUse.set_over(cmapSpec[-1])
    normSpecUse = BoundaryNorm(clevs, cmapSpecUse.N)

    return cmapSpecUse, normSpecUse



def _save_plot(fig: plt.figure, parameter: CoreParameter):
    """Save the plot using the figure object and parameter configs.

    This function creates the output filename to save the plot. It also
    saves each individual subplot if the reference name is an empty string ("").

    Parameters
    ----------
    fig : plt.figure
        The plot figure.
    parameter : CoreParameter
        The CoreParameter with file configurations.
    """
    for f in parameter.output_format:
        f = f.lower().split(".")[-1]
        #fnm = os.path.join(
        #   f'{parameter.spec_type}',
        #    #get_output_dir(parameter.current_set, parameter),
        #    parameter.output_file + "." + f,
        #)
        fnm = f'{parameter.var_id}_{parameter.spec_type}_15N-15N.png'
        plt.savefig(fnm)
        logger.info(f"Plot saved in: {fnm}")

    # Save individual subplots
    if parameter.ref_name == "":
        panels = [PANEL[0]]
    else:
        panels = PANEL

    for f in parameter.output_format_subplot:
        fnm = os.path.join(
            get_output_dir(parameter.current_set, parameter),
            parameter.output_file,
        )
        page = fig.get_size_inches()

        for idx, panel in enumerate(panels):
            # Extent of subplot
            subpage = np.array(panel).reshape(2, 2)
            subpage[1, :] = subpage[0, :] + subpage[1, :]
            subpage = subpage + np.array(BORDER_PADDING).reshape(2, 2)
            subpage = list(((subpage) * page).flatten())  # type: ignore
            extent = Bbox.from_extents(*subpage)

            # Save subplot
            fname = fnm + ".%i." % idx + f
            plt.savefig(fname, bbox_inches=extent)

            orig_fnm = os.path.join(
                './',
                #get_output_dir(parameter.current_set, parameter),
                parameter.output_file,
            )
            fname = orig_fnm + ".%i." % idx + f
            logger.info(f"Sub-plot saved in: {fname}")



def _wave_frequency_plot(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.figure,
    parameter: CoreParameter,
    title: Tuple[str | None, str, str],
    do_zoom: Boolean = False,
):
    """Create wave frequency plot.

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    var : xr.DataArray
        The variable to plot.
    fig : plt.figure
        The figure object to add the subplot to.
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    title : Tuple[str | None, str, str]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    do_zoom: Boolean
    """
    #TODO link var_id
    varName = parameter.var_id
    #varName = 'PRECT'
    PlotDesc = {}
    PlotDesc['spec_raw_sym'] = {"long_name_desc": f"{varName}: log-base10 of lightly smoothed spectral power of component symmetric about equator", "ref_fig_num": "Figure 1"} # Figure number from  Wheeler and Kiladis (1999)
    PlotDesc['spec_raw_asy'] = {"long_name_desc": f"{varName}: log-base10 of lightly smoothed spectral power of component antisymmetric about equator", "ref_fig_num": "Figure 1"}
    PlotDesc['spec_norm_sym'] = {"long_name_desc": f"{varName}: lightly smoothed spectral power of component symmetric about equator, normalized by heavily smoothed background spectrum", "ref_fig_num": "Figure 3"}
    PlotDesc['spec_norm_asy'] = {"long_name_desc": f"{varName}: lightly smoothed spectral power of component antisymmetric about equator, normalized by heavily smoothed background spectrum", "ref_fig_num": "Figure 3"}
    PlotDesc['spec_background'] = {"long_name_desc": f"{varName}: heavily smoothed version of the mean of spectral powers associated with the components symmetric and antisymmetric about equator", "ref_fig_num": "Figure 2"}

    text_offset = 0.005
    fb  = [0, .8]    # frequency bounds for plot
    wnb = [-15, 15]  # zonal wavenumber bounds for plot

    if(max(var['frequency'].values) == 0.5):
        fb = [0, .5]

    if(do_zoom):
        fb  = [0, .18]
        wnb = [-7, 7]

    # get data for dispersion curves:
    equivDepths=[50, 25, 12]
    swfreq,swwn = wf.genDispersionCurves(Ahe=equivDepths)
    # swfreq.shape # -->(6, 3, 50)

    if 'spec_norm' in var.name: # with dispersion curves to plot
        # For n=1 ER waves, allow dispersion curves to touch 0 -- this is for plot aesthetics only
        for i in range(0,3):     # loop 0-->2 for the assumed 3 shallow water dispersion curves for ER waves
          indMinPosFrqER = np.where(swwn[3,i,:] >= 0., swwn[3,i,:], 1e20).argmin()  # index of swwn for least positive wn
          swwn[3,i,indMinPosFrqER],swfreq[3,i,indMinPosFrqER] = 0.,0.    # this sets ER's frequencies to 0. at wavenumber 0.

    swf = np.where(swfreq == 1e20, np.nan, swfreq)
    swk = np.where(swwn == 1e20, np.nan, swwn)


    # Final data refinement:  transpose and trim, set 0 freq to NaN, take log10 for raw, refine metadata
    z  = var.transpose().sel(frequency=slice(*fb), wavenumber=slice(*wnb))
    z.loc[{'frequency':0}] = np.nan

    if 'spec_raw' in var.name:

        east_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(1,3)).sum()
        west_power = z.sel(frequency=slice((1./96.),(1./24.)), wavenumber=slice(-3,-1)).sum()
        ew_ratio   = east_power / west_power
        print("\neast_power: %12.5f" % east_power)
        print("west_power: %12.5f" % west_power)
        print("ew_ratio: %12.5f\n" % ew_ratio)

        z = np.log10(z)

    if 'spec_background' in var.name:
          z = np.log10(z)

    z.attrs["long_name"] = PlotDesc[var.name]["long_name_desc"]
    z.attrs["method"] = f"Follows {PlotDesc[var.name]['ref_fig_num']} methods of Wheeler and Kiladis (1999; https://doi.org/10.1175/1520-0469(1999)056<0374:CCEWAO>2.0.CO;2)"

    if 'spec_raw' in var.name:

        z.attrs["ew_ratio_method"] = "Sum of raw (not log10) symmetric spectral power for ZWNs +/- 1-3, periods 24-96 days"
        z.attrs["east_power"] = east_power.values
        z.attrs["west_power"] = west_power.values
        z.attrs["ew_ratio"]   = ew_ratio.values

#    # TODO Save plotted data z to file as xArray data array
#    dataDesc = f"spectral_power_{srcID}_{vari}_{spec_type}_{component}_200101_201412"
#    z.to_netcdf(outDataDir + "/"+ dataDesc + ".nc")

    #fig, ax = plt.subplots()
    ax = fig.add_axes(PANEL[subplot_num])
 
    kmesh0, vmesh0 = np.meshgrid(z['wavenumber'], z['frequency'])
    #img = ax.contourf(kmesh0, vmesh0, z, levels=np.linspace(0.2, 3.0, 16), cmap='Spectral_r',  extend='both')

    # for test and ref:
    if subplot_num < 2:
        if 'spec_norm' in var.name:
            contour_level_spec = CONTOUR_LEVS_SPEC_NORM
            cmapSpecUse, normSpecUse = create_colormap_clevs(CMAP_SPEC_NORM, CONTOUR_LEVS_SPEC_NORM)
        else:
            contour_level_spec = CONTOUR_LEVS_SPEC_RAW
            cmapSpecUse, normSpecUse = create_colormap_clevs(CMAP_SPEC_RAW, CONTOUR_LEVS_SPEC_RAW)
        img = ax.contourf(kmesh0, vmesh0, z, levels=contour_level_spec, cmap=cmapSpecUse,  norm=normSpecUse, extend='both')
        img2 = ax.contour(kmesh0, vmesh0, z, levels=contour_level_spec, linewidths=1., linestyles='solid', colors='gray', alpha=0.7)
    
    # for diff ratio
    if subplot_num == 2:
        # TODO refine color bar
        if 'spec_norm' in var.name:
            contour_level_spec = CONTOUR_LEVS_SPEC_NORM_DIFF
        else:
            contour_level_spec = CONTOUR_LEVS_SPEC_RAW_DIFF
        cmapSpecUse = 'seismic'

        img = ax.contourf(kmesh0, vmesh0, z, levels=contour_level_spec, cmap=cmapSpecUse, extend='both')
        img2 = ax.contour(kmesh0, vmesh0, z, levels=contour_level_spec, linewidths=1., linestyles='solid', colors='gray', alpha=0.7)

    ax.axvline(0, linestyle='dashed', color='dimgray', linewidth=1.0, alpha=0.60)
    if( (1./30.) < fb[1] ):
      ax.axhline((1./30.), linestyle='dashed', color='dimgray', alpha=0.80)
      ax.text(wnb[0]+1,(1./30.)+text_offset,'30 days',color='dimgray', alpha=0.80)
    if( (1./6.) < fb[1] ):
      ax.axhline((1./6.), linestyle='dashed', color='dimgray', alpha=0.80)
      ax.text(wnb[0]+1,(1./6.)+text_offset,'6 days',color='dimgray', alpha=0.80)
    if( (1./3.) < fb[1] ):
      ax.axhline((1./3.), linestyle='dashed', color='dimgray', alpha=0.80)
      ax.text(wnb[0]+1,(1./3.)+text_offset,'3 days',color='dimgray', alpha=0.80)
    for ii in range(3,6):
        ax.plot(swk[ii, 0,:], swf[ii,0,:], color='black', linewidth=1.5, alpha=0.80)
        ax.plot(swk[ii, 1,:], swf[ii,1,:], color='black', linewidth=1.5, alpha=0.80)
        ax.plot(swk[ii, 2,:], swf[ii,2,:], color='black', linewidth=1.5, alpha=0.80)
    ax.set_xlim(wnb)
    ax.set_ylim(fb)
    #ax.set_title(varName + ": Log ${\sum_{15^{\circ}S}^{15^{\circ}N} Power_{SYM}}$")   # Version w/ LaTeX
    if 'spec_raw' in var.name:
        ax.set_title(f"{varName}: Log{{Sum(Power) from 15°S-15°N}}\n")       # Version w/o LaTeX
    elif 'spec_norm' in var.name:
        ax.set_title(f"{varName}: {{Sum(Power) from 15°S-15°N}}/Background\n")          # Version w/o LaTeX
    else:
        ax.set_title(f"{varName}: Log{{Smoothed Background Power}}\n")

    ax.set_title('model', loc='left')
    print('*****',var)
    ax.set_title(f"{var.component}", loc='right')

    if 'spec_norm' in var.name:
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
    fig.text(
	    PANEL[subplot_num][0] - 0.03 ,
            PANEL[subplot_num][1] - 0.03,  "Westward", fontsize=11)
    fig.text(
	    PANEL[subplot_num][0] + 0.35 ,
            PANEL[subplot_num][1] -0.03 ,  "Eastward", fontsize=11)
    fig.colorbar(img)
#    # Save fig
#    fig.savefig(outDataDir + "/"+ dataDesc + "_plot.png", bbox_inches='tight', dpi=300)
#
#    print("Plot file created: %s\n" % outDataDir + "/"+ dataDesc + "_plot.png")




def plot(
    parameter: CoreParameter,
    da_test: xr.DataArray,
    da_ref: xr.DataArray | None,
    da_diff: xr.DataArray | None,
    do_zoom: Boolean = False,
):
    """Plot the variable's metrics generated for the lat_lon set.

    Parameters
    ----------
    parameter : CoreParameter
        The CoreParameter object containing plot configurations.
    da_test : xr.DataArray
        The test data.
    da_ref : xr.DataArray | None
        The optional reference data.
    ds_diff : xr.DataArray | None
        The difference between ``ds_test`` and ``ds_ref``.
    """
    #fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig = plt.figure(figsize=[8.5, 12.0], dpi=300)
    #fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)


    _wave_frequency_plot(
        0,
        da_test,
        fig,
        parameter,
        title=(parameter.test_name_yrs, parameter.test_title),  # type: ignore
        do_zoom=do_zoom,
    )


    _wave_frequency_plot(
        1,
        da_ref,
        fig,
        parameter,
        title=(parameter.ref_name_yrs, parameter.reference_title),  # type: ignore
        do_zoom=do_zoom,
    )


    _wave_frequency_plot(
        2,
        da_diff,
        fig,
        parameter,
        title=(None, parameter.diff_title),  # type: ignore
        do_zoom=do_zoom,
    )

#TODO: save plot:NameError: name 'get_output_dir' is not defined
    _save_plot(fig, parameter)

    plt.close()
