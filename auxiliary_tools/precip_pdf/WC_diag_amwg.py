#!/usr/local/uvcdat/1.3.1/bin/python

# ======================================================

# Library with functions for diagnostics of the water cycle
# Created by Chris Terai (terai1@llnl.gov)
# Before running: load e3sm_unified package
# on cori it is: source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified.sh

# ======================================================


#from unidata import udunits
import cdutil
import cdutil.times, numpy
from numbers import Number
from pprint import pprint
import cdms2
import genutil
import MV2
#many of the following are for writing attributes to netcdf
import cdat_info,cdtime,code,datetime,gc,inspect,os,re,string,sys,pytz
from socket import gethostname
#from string import replace
# load potentially relevant libraries

import matplotlib.pyplot as plt
import numpy.ma as ma
import matplotlib.colors as colors
import cartopy.crs as ccrs
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter




def create_precip_PDF_netcdf(mv, model_handle, model_filename_str, path='',model_time='0001'):
    """ If a netcdf file containing the probability distribution function (PDF) 
    does not exist for the model data, this script takes the data and 
    creates a netcdf with the frequency and amount PDFs based on the 
    specified data. Name of the file will follow other climo files with the format: 
    filename = case + variable + "_PDF_climo.nc" 
    'path' defines the location on which to write the files
    
    Currently, the bin edge values and the scale of the bin widths are specified in the script.
    This feature can and may be changed in the future.
    """
    
    #case= model_handle.attributes['case']
    #simulation_name,grid_resolution,version_name,others_name=model_filename_str.split('_')
    #ensemble_name,model_name,timetype_name,model_time,file_type=others_name.split('.')

    #model_time='0005-01-01'
    variablename=mv.id
    vid_m=''.join([variablename,'FREQPDF'])
    vid2_m=''.join([variablename,'AMNTPDF'])
    vid3_m='bincenter'
    
    print("calculating pdf for ",variablename)
    #specify bin edges of PDF
    #specify minimum precipitation threshold
    threshold=0.1     # in mm d-1
    bin_increase=1.07 # multiplicative factor
    binedges_input=threshold * bin_increase ** (numpy.arange(0,130))   #these are the bin edges
    [bincenter, output_mapped_freqpdf, output_mapped_amntpdf] = create_amount_freq_PDF(mv,binedges_input, binwidthtype='logarithmic', bincentertype='arithmetic', vid=vid_m, vid2=vid2_m,vid3=vid3_m)
    
    print("writing pdf file for ",model_time,variablename)
     
    outputfilename=model_filename_str
    #outputfilename='_'.join([case,model_time,variablename,'PDF.nc'])
    
    
    f_out = cdms2.open(os.path.join(path,outputfilename),'w') 
    att_keys = model_handle.attributes.keys()
    att_dic = {}
    #for i in range(len(att_keys)):
    #    att_dic[i]=att_keys[i],model_handle.attributes[att_keys[i]]
    #    to_out = att_dic[i]
    #    if not hasattr(f_out,to_out[0]):
    #        setattr(f_out,to_out[0],to_out[1])

    #parts taken from globAttWrite by Paul Durack - write out who, what, where when, the data was processed
    #   Comment out the following if pytz is not installed on machine
    local                       = pytz.timezone("America/Los_Angeles")
    time_now                    = datetime.datetime.now();
    #local_time_now              = time_now.replace(tzinfo = local)
    #utc_time_now                = local_time_now.astimezone(pytz.utc)
    #time_format                 = utc_time_now.strftime("%d-%b-%Y %H:%M:%S %p")
    #f_out.institution     = "Program for Climate Model Diagnosis and Intercomparison (LLNL)"
    f_out.data_contact    = "Chris Terai; terai1@llnl.gov; +1 925 423 7972"
    #f_out.history         = "".join(['File processed: ',time_format,' UTC; San Francisco, CA, USA ',f_out.history])
    f_out.host            = "".join([gethostname(),'; CDAT version: ',".".join(["%s" % el for el in cdat_info.version()])])
    f_out.write(bincenter) ;
    f_out.write(output_mapped_freqpdf) ;
    f_out.write(output_mapped_amntpdf) ;
    f_out.close()
    print("pdf file written for ",outputfilename," ",variablename)
    print("  **************    ")

def create_precip_PDF_netcdf_ConvTest(mv, model_handle, model_prefix, path=''):
    """ If a netcdf file containing the probability distribution function (PDF) 
    does not exist for the model data, this script takes the data and 
    creates a netcdf with the frequency and amount PDFs based on the 
    specified data. Name of the file will follow other climo files with the format: 
    filename = case + variable + "_PDF_climo.nc" 
    'path' defines the location on which to write the files
    
    Currently, the bin edge values and the scale of the bin widths are specified in the script.
    This feature can and may be changed in the future.
    """
    
    variablename=mv.id
    vid_m=''.join([variablename,'FREQPDF'])
    vid2_m=''.join([variablename,'AMNTPDF'])
    vid3_m='bincenter'
    
    print("calculating pdf for ",model_prefix,variablename)
    #specify bin edges of PDF
    #specify minimum precipitation threshold
    threshold=0.1     # in mm d-1
    bin_increase=1.07 # multiplicative factor
    binedges_input=threshold * bin_increase ** (numpy.arange(0,130))   #these are the bin edges
    [bincenter, output_mapped_freqpdf, output_mapped_amntpdf] = create_amount_freq_PDF(mv,binedges_input, binwidthtype='logarithmic', bincentertype='arithmetic', vid=vid_m, vid2=vid2_m,vid3=vid3_m)
    
    print("writing pdf file for ",model_prefix,variablename)
     
    
    outputfilename='_'.join([model_prefix,variablename,'PDF.nc'])
    
    
    f_out = cdms2.open(os.path.join(path,outputfilename),'w') 
    att_keys = model_handle.attributes.keys()
    att_dic = {}
    for i in range(len(att_keys)):
        att_dic[i]=att_keys[i],model_handle.attributes[att_keys[i]]
        to_out = att_dic[i]
        if not hasattr(f_out,to_out[0]):
            setattr(f_out,to_out[0],to_out[1])

    #parts taken from globAttWrite by Paul Durack - write out who, what, where when, the data was processed
    #   Comment out the following if pytz is not installed on machine
    local                       = pytz.timezone("America/Los_Angeles")
    time_now                    = datetime.datetime.now();
    #local_time_now              = time_now.replace(tzinfo = local)
    #utc_time_now                = local_time_now.astimezone(pytz.utc)
    #time_format                 = utc_time_now.strftime("%d-%b-%Y %H:%M:%S %p")
    #f_out.institution     = "Program for Climate Model Diagnosis and Intercomparison (LLNL)"
    f_out.data_contact    = "Chris Terai; terai1@llnl.gov; +1 925 422 8830"
    #f_out.history         = "".join(['File processed: ',time_format,' UTC; San Francisco, CA, USA '])
    f_out.host            = "".join([gethostname(),'; CDAT version: ',".".join(["%s" % el for el in cdat_info.version()])])
    f_out.write(bincenter) ;
    f_out.write(output_mapped_freqpdf) ;
    f_out.write(output_mapped_amntpdf) ;
    f_out.close()
    print("pdf file written for ",model_prefix,variablename)
    print("  **************    ")

def create_amount_freq_PDF(mv, binedges, binwidthtype=None, bincentertype=None, vid=None, vid2=None, vid3=None):
    """Takes in geospatial data (mv) with dimensions of lat, lon, and time, and 
    creates a PDF of the mv based on binedges at each lat/lon grid point. 
    binedges defines the edges of the bin, except for the bin with maximum value,
    where it is open. 
    
    binwidth option allows user to define whether the PDF is scaled by the 
    'arithmetic' bin width or the 'logarithmic' bin width. (dN/dx vs. dN/dlogx)
    default: 'logarithmic'
    
    The bincenter option allows one to use the 'geometric' 
    or the 'arithmetic' mean to define the edge. The bin with maximum value will have
    a bin center equidistant from the max bin edge as from the center of the previous bin. 
    
    PDFs will not be normalized over the histogram, but over all available data.
    For example, if there is data below the minimum bin edge, then the PDF will be 
    normalized, having included the data that lies below the minimum bin edge.
    
    vid = variable ID, which will typically be the variable name of the when output as a netcdf
    """
    
    
    #Step 1 input data and figure out the time dimensions in the data
    if vid is None:
        vid = mv.id
        vid2 = ''.join([mv.id,'2'])
    
    
        
    #Do get domain and find which axis corresponds to time, lat and lon
    time_index=mv.getAxisIndex('time')
    lat_index=mv.getAxisIndex('lat')
    lon_index=mv.getAxisIndex('lon')
    #obtain long_name, standard_name, typecode of the variable to eventually feed into output variables
    var_long_name=mv.long_name
    mv_typecode=mv.typecode()
    mv_lat=mv.getAxis(lat_index)
    mv_lon=mv.getAxis(lon_index)
    mv_att=mv.attributes
    mv_grid=mv.getGrid()
    
    
    #Step 2 loop over the bin widths and add up the number of data points in each bin
    
    #Create an array with the shape of (lat,lon,binedges,and corresponding bincenter)
    mapped_precip_freqpdf=numpy.zeros((mv.shape[lat_index],mv.shape[lon_index],len(binedges)))
    mapped_precip_amntpdf=numpy.zeros((mv.shape[lat_index],mv.shape[lon_index],len(binedges)))
    bincenter=numpy.zeros(len(binedges))
    
    #Count up total first
    counts_index=MV2.greater_equal(mv,0.)
    data_counts=numpy.zeros((mv.shape))
    data_counts[counts_index]=1.
    counts_total=numpy.sum(data_counts,axis=time_index)
    
    #specify what the binmean, bincenter, and binwidths are based on log and arith scaling
    binwidth=numpy.zeros(len(binedges))
    bincenter=numpy.zeros(len(binedges))
    binmean=numpy.zeros(len(binedges))
    
    #Calculate bin mean for amount PDF
    binmean[:-1]=(binedges[1:]+binedges[:-1])/2
    binmean[-1]=binedges[-1] + (binedges[-1]-binedges[-2])/2
    #Calculate bin width based on type
    if binwidthtype is 'arithmetic':
        binwidth[:-1]=binedges[1:]-binedges[:-1]
        binwidth[-1]=binedges[-1]-binedges[-2]
    elif binwidthtype is 'logarithmic' or binwidthtype is None:
        binwidth[:-1]=numpy.log10(binedges[1:]/binedges[:-1])
        binwidth[-1]=numpy.log10(binedges[-1]/binedges[-2])
    #Calculate bin center based on type
    if bincentertype is 'arithmetic' or bincentertype is None:
        bincenter=binmean
    elif bincentertype is 'geometric':
        bincenter[:-1]=numpy.sqrt(binedges[1:]*binedges[:-1])
        bincenter[-1]=binedges[-1] + (binedges[-1]-binedges[-2])/2
    
    #Count up the number of days of precip in each precip bin **Most work done here
    for i in range(len(binedges)):
        precip_index=numpy.ones(mv.shape)  #locate the index where precip rate is between the bin edges
        toolow_index=mv<binedges[i]        
        precip_index[toolow_index]=0
        if i!=(len(binedges)-1):
            toohigh_index=mv>=binedges[i+1]
            precip_index[toohigh_index]=0
        precip_total=numpy.sum(precip_index,axis=time_index)
        precip_fraction=numpy.divide(precip_total,counts_total)
        
        precip_freqpdf=precip_fraction/binwidth[i]
        precip_amntpdf=precip_fraction/binwidth[i]*binmean[i]
        mapped_precip_freqpdf[:,:,i]=precip_freqpdf
        mapped_precip_amntpdf[:,:,i]=precip_amntpdf
        precip_freqpdf=None
        precip_amntpdf=None
        
    
    
    #Step 3 attach all necessary attributes to data (create data as a transient variable)
    #First, specify a new axis for the PDF
    binbound_all=numpy.append(binedges,numpy.max(mv))
    binbounds=numpy.zeros((len(binedges),2))
    binbounds[:,0]=binbound_all[:-1]
    binbounds[:,1]=binbound_all[1:]
    mv_hist=cdms2.createAxis(bincenter,bounds=binbounds,id='binvalue') #mv_hist is the axes for precip rate
    ouput_mapped_precip_freqpdf=cdms2.createVariable(mapped_precip_freqpdf,typecode=mv_typecode,
                                             grid=mv_grid,axes=[mv_lat,mv_lon,mv_hist],attributes=mv_att,id=vid)
    
    ouput_mapped_precip_amntpdf=cdms2.createVariable(mapped_precip_amntpdf,typecode=mv_typecode,
                                             grid=mv_grid,axes=[mv_lat,mv_lon,mv_hist],attributes=mv_att,id=vid2)
    bincenter=cdms2.createVariable(bincenter,typecode=mv_typecode, axes=[mv_hist],attributes=mv_att,id=vid3)
                                       
    ouput_mapped_precip_freqpdf.units='frequency'
    ouput_mapped_precip_amntpdf.units='amount'
    ouput_mapped_precip_freqpdf.long_name=''.join(['Frequency as a function of ',var_long_name])
    ouput_mapped_precip_amntpdf.long_name=''.join(['Amount as a function of ',var_long_name])
    try:
        var_standard_name=mv.standard_name
        ouput_mapped_precip_freqpdf.standard_name=''.join([var_standard_name,'_frequency'])
        ouput_mapped_precip_amntpdf.standard_name=''.join([var_standard_name,'_amount'])
    except:
        print("Variable does not have standard_name")
    #Step 4 output data
    return bincenter, ouput_mapped_precip_freqpdf, ouput_mapped_precip_amntpdf

# Set of definitions needed to run the plotting                                                                                                                                     
def cross_section(var,lat,lon,lat_lim,lon_lim,lat_offset,direction):
    # Will take in variable over location with lat and lon                                                                                                                          
    # limits and creates cross_section with equal latitude                                                                                                                          
    # steps of lat_offset*2                                                                                                                                                         
    # assumes that there's some width to the latitude (lat_lim[1]>lat_lim[0])                                                                                                       

    # assumes var has shape of [lev,lat,lon]                                                                                                                                        
    # Figure out the dimensions of the variable                                                                                                                                     
    lev_coord_length=var.shape[0]
    lat_coord_length=var.shape[1]
    lon_coord_length=var.shape[2]

    #count number of columns in transect (cross-section)                                                                                                                            
    lat_slice_count=int(numpy.ceil(numpy.absolute(lat_lim[1]-lat_lim[0])/(lat_offset*2)))
    print('lat slice count is ',lat_slice_count)

    lat_cent=numpy.zeros((lat_slice_count-1,1))
    lon_cent=numpy.zeros((lat_slice_count-1,1))
    #determine lat_center values of transect                                                                                                                                        
    lat_cent=lat_lim[0]+lat_offset+numpy.arange(lat_slice_count-1)*lat_offset*2
    lon_offset=((lon_lim[1]-lon_lim[0])/lat_slice_count)/2
    print('lon offset is ',lon_offset)
    #determine the lon center values of transect                                                                                                                                    
    lon_cent=lon_lim[0]+lon_offset+numpy.arange(lat_slice_count-1)*lon_offset*2

    lat_edge=numpy.zeros((lat_slice_count-1,1))
    lon_edge=numpy.zeros((lat_slice_count-1,1))

    lat_edge=lat_lim[0]+numpy.arange(lat_slice_count)*lat_offset*2
    lon_edge=lon_lim[0]+numpy.arange(lat_slice_count)*lon_offset*2
    lon_width=2+numpy.abs(lon_offset)
    print('lon width is ',lon_width)


    var_crosssection=numpy.zeros((lev_coord_length,lat_slice_count-1))
    var_crosssection_pcolor=numpy.zeros((lev_coord_length,lat_slice_count))
    for i in numpy.arange(lat_slice_count-1):
        lat_arg=numpy.squeeze(numpy.argwhere((lat>lat_cent[i]-lat_offset) & (lat<lat_cent[i]+lat_offset)))
        lon_arg=numpy.squeeze(numpy.argwhere((lon<lon_cent[i]+lon_width) & (lon>lon_cent[i]-lon_width)))
        var_lonslice=var[:,:,lon_arg]
        var_latslice=var_lonslice[:,lat_arg,:]
        if (len(var_latslice.shape)>2):
            var_column=numpy.mean(var_latslice,axis=(1,2))
        else:
            var_column=numpy.mean(var_latslice,axis=(1))
        var_crosssection[:,i]=var_column
        var_crosssection_pcolor[:,i]=var_column
    return var_crosssection,var_crosssection_pcolor,lat_cent,lon_cent,lat_edge,lon_edge

# Set of definitions needed to run the plotting
from matplotlib.colors import LinearSegmentedColormap

def plot_panel(n, fig, proj, var, clevels, cmap,
               title, stats=None):
    panel = [(0.1691, 0.6810, 0.6465, 0.2258),
         (0.1691, 0.3961, 0.6465, 0.2258),
         (0.1691, 0.1112, 0.6465, 0.2258),
         ]
    plotTitle = {'fontsize': 11.5}
    plotSideTitle = {'fontsize': 9.5}
    #var = add_cyclic(var)
    lon = var.getLongitude()
    lat = var.getLatitude()
    var = ma.squeeze(var.asma())

    # Contour levels
    levels = None
    norm = None
    if len(clevels) > 0:
        levels = [-1.0e8] + clevels + [1.0e8]
        norm = colors.BoundaryNorm(boundaries=levels, ncolors=256)

    #ax.set_global()
    region_str = 'global'
    #region = regions_specs[region_str]
    global_domain = True
    full_lon = True
    #if 'domain' in region.keys():
    #    # Get domain to plot
    #    domain = region['domain']
    #    global_domain = False
    #else:
        # Assume global domain
    domain = cdutil.region.domain(latitude=(-90., 90, 'ccb'))
    kargs = domain.components()[0].kargs
    lon_west, lon_east, lat_south, lat_north = (0, 360, -90, 90)
    if 'longitude' in kargs:
        full_lon = False
        lon_west, lon_east, _ = kargs['longitude']
        # Note cartopy Problem with gridlines across the dateline:https://github.com/SciTools/cartopy/issues/821. Region cross dateline is not supported yet.
        if lon_west>180 and lon_east>180:
            lon_west = lon_west - 360
            lon_east = lon_east - 360
            
    if 'latitude' in kargs:
        lat_south, lat_north, _ = kargs['latitude']
    lon_covered = lon_east - lon_west
    lon_step = determine_tick_step(lon_covered)
    xticks = numpy.arange(lon_west, lon_east, lon_step)
    # Subtract 0.50 to get 0 W to show up on the right side of the plot.
    # If less than 0.50 is subtracted, then 0 W will overlap 0 E on the left side of the plot.
    # If a number is added, then the value won't show up at all.
    if global_domain or full_lon:
        xticks = numpy.append(xticks, lon_east-0.50)
        proj = ccrs.PlateCarree(central_longitude=180)
    else:
        xticks = numpy.append(xticks, lon_east)
    lat_covered = lat_north - lat_south
    lat_step = determine_tick_step(lat_covered)
    yticks = numpy.arange(lat_south, lat_north, lat_step)
    yticks = numpy.append(yticks, lat_north)

    # Contour plot
    ax = fig.add_axes(panel[n], projection=proj)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=proj)
    cmap=plt.get_cmap(cmap)
    #cmap = get_colormap(cmap)
    p1 = ax.contourf(lon, lat, var,
                     transform=ccrs.PlateCarree(),
                     norm=norm,
                     levels=levels,
                     cmap=cmap,
                     extend='both',
                     )
    
    #ax.set_aspect('auto')
    # Full world would be aspect 360/(2*180) = 1
    ax.set_aspect((lon_east - lon_west)/(2*(lat_north - lat_south)))
    ax.coastlines(lw=0.3)
    if not global_domain and 'RRM' in region_str:
        ax.coastlines(resolution='50m', color='black', linewidth=1)
        state_borders = cfeature.NaturalEarthFeature(category='cultural', name='admin_1_states_provinces_lakes', scale='50m', facecolor='none')
        ax.add_feature(state_borders, edgecolor='black')
    if title[0] is not None:
        ax.set_title(title[0], loc='left', fontdict=plotSideTitle)
    if title[1] is not None:
        ax.set_title(title[1], fontdict=plotTitle)
    if title[2] is not None:
        ax.set_title(title[2], loc='right', fontdict=plotSideTitle)
    ax.set_xticks(xticks, crs=ccrs.PlateCarree())
    ax.set_yticks(yticks, crs=ccrs.PlateCarree())
    lon_formatter = LongitudeFormatter(
        zero_direction_label=True, number_format='.0f')
    lat_formatter = LatitudeFormatter()
    ax.xaxis.set_major_formatter(lon_formatter)
    ax.yaxis.set_major_formatter(lat_formatter)
    ax.tick_params(labelsize=8.0, direction='out', width=1)
    ax.xaxis.set_ticks_position('bottom')
    ax.yaxis.set_ticks_position('left')

    # Color bar
    cbax = fig.add_axes(
        (panel[n][0] + 0.6635, panel[n][1] + 0.0215, 0.0326, 0.1792))
    cbar = fig.colorbar(p1, cax=cbax)
    w, h = get_ax_size(fig, cbax)

    if levels is None:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    else:
        maxval = numpy.amax(numpy.absolute(levels[1:-1]))
        if maxval < 10.0:
            fmt = "%5.2f"
            pad = 25
        elif maxval < 100.0:
            fmt = "%5.1f"
            pad = 25
        else:
            fmt = "%6.1f"
            pad = 30
        cbar.set_ticks(levels[1:-1])
        labels = [fmt % l for l in levels[1:-1]]
        cbar.ax.set_yticklabels(labels, ha='right')
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)

def plot(reference, test, diff, figsize_in, contour_levels, test_colormap, test_name_yrs, \
         test_title, units, ref_name_yrs, reference_title, diff_levels, diff_colormap, \
         diff_title, main_title, plt_fig_name):

    # Create figure, projection
    fig = plt.figure(figsize=figsize_in, dpi=300)
    proj = ccrs.PlateCarree()
    # First two panels
    #min1 = metrics_dict['test']['min']
    #mean1 = metrics_dict['test']['mean']
    #max1 = metrics_dict['test']['max']

    plot_panel(0, fig, proj, test, contour_levels, test_colormap,
               (test_name_yrs, test_title, units))#, stats=(max1, mean1, min1))

    #min2 = metrics_dict['ref']['min']
    #mean2 = metrics_dict['ref']['mean']
    #max2 = metrics_dict['ref']['max']
    plot_panel(1, fig, proj, reference, contour_levels, test_colormap,
               (ref_name_yrs, reference_title, units))#, stats=(max2, mean2, min2))

    # Third panel
    #min3 = metrics_dict['diff']['min']
    #mean3 = metrics_dict['diff']['mean']
    #max3 = metrics_dict['diff']['max']
    #r = metrics_dict['misc']['rmse']
    #c = metrics_dict['misc']['corr']
    plot_panel(2, fig, proj, diff, diff_levels, diff_colormap,
               (None, diff_title, None))#, stats=(max3, mean3, min3, r, c))

    # Figure title
    fig.suptitle(main_title, x=0.5, y=0.96, fontsize=18)
    # Save figure    
    plt.savefig(plt_fig_name)
    plt.close()


def determine_tick_step(degrees_covered):
    if degrees_covered > 180:
        return 60
    if degrees_covered > 60:
        return 30
    elif degrees_covered > 30:
        return 10
    elif degrees_covered > 20:
        return 5
    else:
        return 1

def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0], lon[0] + 360.0, 'coe'))

def get_ax_size(fig, ax):
    bbox = ax.get_window_extent().transformed(fig.dpi_scale_trans.inverted())
    width, height = bbox.width, bbox.height
    width *= fig.dpi
    height *= fig.dpi
    return width, height

# Plotting of the precip PDF
#   Global and tropical mean
def Global_PDF_values(file_handle,variable):
    freqpdfvarname=''.join([variable,'FREQPDF'])
    amntpdfvarname=''.join([variable,'AMNTPDF'])
    PRECT_FREQPDF=file_handle(freqpdfvarname,lat=[-90,90],binvalue=[0.1,600.])
    PRECT_AMNTPDF=file_handle(amntpdfvarname,lat=[-90,90],binvalue=[0.1,600.])
    bincenter=file_handle('bincenter',binvalue=[0.1,600.])
    global_mean_FREQPDF=cdutil.averager(PRECT_FREQPDF,axis='xy')
    global_mean_AMNTPDF=cdutil.averager(PRECT_AMNTPDF,axis='xy')
    bincenter=numpy.array(bincenter)
    globalmeanFREQPDF=numpy.array(global_mean_FREQPDF)
    globalmeanAMNTPDF=numpy.array(global_mean_AMNTPDF)
    return bincenter,globalmeanFREQPDF,globalmeanAMNTPDF

# For specific region, change the lat bounds below from -30,30 and add lon=[XX,YY], 
def Tropical_PDF_values(file_handle,variable):
    freqpdfvarname=''.join([variable,'FREQPDF'])
    amntpdfvarname=''.join([variable,'AMNTPDF'])
    PRECT_FREQPDF=file_handle(freqpdfvarname,lat=[-30,30],binvalue=[0.1,600.])
    PRECT_AMNTPDF=file_handle(amntpdfvarname,lat=[-30,30],binvalue=[0.1,600.])
    bincenter=file_handle('bincenter',binvalue=[0.1,600.])
    global_mean_FREQPDF=cdutil.averager(PRECT_FREQPDF,axis='xy')
    global_mean_AMNTPDF=cdutil.averager(PRECT_AMNTPDF,axis='xy')
    bincenter=numpy.array(bincenter)
    globalmeanFREQPDF=numpy.array(global_mean_FREQPDF)
    globalmeanAMNTPDF=numpy.array(global_mean_AMNTPDF)
    return bincenter,globalmeanFREQPDF,globalmeanAMNTPDF

def CONUS_PDF_values(file_handle,variable):
    freqpdfvarname=''.join([variable,'FREQPDF'])
    amntpdfvarname=''.join([variable,'AMNTPDF'])
    PRECT_FREQPDF=file_handle(freqpdfvarname,lat=[35,49],lon=[-125,-75.],binvalue=[0.1,600.])
    PRECT_AMNTPDF=file_handle(amntpdfvarname,lat=[35,49],lon=[-125,-75.],binvalue=[0.1,600.])
    bincenter=file_handle('bincenter',binvalue=[0.1,600.])
    global_mean_FREQPDF=cdutil.averager(PRECT_FREQPDF,axis='xy')
    global_mean_AMNTPDF=cdutil.averager(PRECT_AMNTPDF,axis='xy')
    bincenter=numpy.array(bincenter)
    globalmeanFREQPDF=numpy.array(global_mean_FREQPDF)
    globalmeanAMNTPDF=numpy.array(global_mean_AMNTPDF)
    return bincenter,globalmeanFREQPDF,globalmeanAMNTPDF
