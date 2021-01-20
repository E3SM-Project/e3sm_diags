import os
from .utils import add_header, h1_to_h3
from .default_viewer import create_metadata
from cdp.cdp_viewer import OutputViewer
import csv
#from varid_dict import varid_longname
from collections import OrderedDict

varid_longname = {'tas':'Surface Temperature (C)','pr':'Precipitation (mm/day)','clt':'Total Cloud Fraction (%)','hurs':'Rel. Humidity (%)','hfss':'Sensible Heat Flux (W/m2)','hfls':'Latent Heat Flux(W/m2)','rlus':'Upwelling LW (W/m2)','rlds':'Downwelling LW (W/m2)','rsus':'Upwelling SW (W/m2)','rsds':'Downwelling SW (W/m2)','ps':'Surface Pressure (Pa)','prw':'Preciptable Water (mm)','cllvi':'Liquid Water Path (mm)','albedo':'Surface Albedo','cl_p':'Cloud Fraction (%)','mrsos':'top 10 cm soil moisture content (mm)','od550aer':'aerosol optical depth'}


def create_viewer(root_dir, parameters):
    """
    Given a set of parameter for a certain set of diagnostics,
    create a single page.

    Return the title and url for this page.
    """
    viewer = OutputViewer(path=root_dir)
    # The name that's displayed on the viewer.
    set_name = 'arm_diags'
    display_name = 'Diagnostics at ARM stations'
    # The title of the colums on the webpage.
    # Appears in the second and third columns of the bolded rows.
    cols = ['Description', 'Plot']
    viewer.add_page(display_name, short_name=set_name, columns=cols)
    viewer.add_group('Annual Cycle')
    viewer.add_row('Variables at SGP')
    viewer.add_col('Annual cycles')
    viewer.add_col('Line plots')
    #viewer.add_row('Variables at Darwin')
    #viewer.add_col('Annual cycles')
    #viewer.add_col('Line plots')

    viewer.add_group('Seasonal Mean')
    viewer.add_row('All variables at SGP')
    viewer.add_col('Seasonal Mean')
    viewer.add_col('ANN')
    viewer.add_col('DJF')
    viewer.add_col('JJA')
    viewer.add_row('All variables at Darwin')
    viewer.add_col('Seasonal Mean')
    viewer.add_col('ANN')
    viewer.add_col('DJF')
    viewer.add_col('JJA')
    viewer.add_col('Annual cycles')

    viewer.add_group('Diurnal Cycle')
    viewer.add_row('PRECT at SGP') 
    viewer.add_col('Diurnal cycle of precipitation')
    viewer.add_col('plot')

    viewer.add_group('Diurnal Cycle/Annual Cycle')
    viewer.add_row('CLOUD at SGP') 
    viewer.add_col('Diurnal cycle and Annual cycle of CLOUD')
    viewer.add_col('plot')

    url = viewer.generate_page()
    add_header(root_dir, os.path.join(root_dir, url), parameters)
    h1_to_h3(os.path.join(root_dir, url))

    return display_name, url

    
#    viewer.add_group('Results of addition')
#    viewer.add_row('Result of 1 + 1')
#    viewer.add_col('Some description for add')
#    viewer.add_col('output.png', is_file=True)
#    viewer.add_row('Another Result')
#    viewer.add_col('Another description for add')
#    viewer.add_col('output.png', is_file=True)
#    
#    viewer.add_group('Results of subtraction')
#    viewer.add_row('Some Result')
#    viewer.add_col('Some description for sub')
#    viewer.add_col('output.png', is_file=True)
#    
#    viewer.generate_viewer()
#    #url = viewer.generate_page()
#    #add_header(root_dir, os.path.join(root_dir, url), parameters)
#    #h1_to_h3(os.path.join(root_dir, url))
#
#    return display_name, url



#def create_viewer(root_dir,parameter):
#    """Creat the main html page hosting all sets of diagnostics"""
#    #output_path = parameter.output_path
#    test_name = parameter[0].short_test_name if parameter[0].short_test_name else parameter[0].test_name
#    output_path = root_dir
#    #set_name = 'arm_diags_annual_cycle'
#    #display_name = 'ARM Diags Annual Cycle'
#    new_dir = os.path.join(root_dir,set_name)
#    if not os.path.exists(new_dir):
#        os.makedirs(new_dir)
#    index_file = os.path.join(new_dir, 'index.html')
#    url = os.path.join(set_name,'index.html')
#
#    f = open(index_file,'w')
#    
#    #annual_cycle_html(parameter)
#
#   
#    message = """<html>
#    <head>
#    <TITLE>ARM Diagnostics Plots</TITLE>
#    </head>
#    </table>
#    </b></font>
#    <p>
#    <b>ARM Metrics and Diagnostics Package</b>
#    <p>
#    <b>Model: """+test_name+"""</b>
#    <hr noshade size=2 size="100%">
#    <TABLE width='1550' >
#    <TR>
#    <TD>
#      <TH ALIGN=left VALIGN=top>
#      <font color=blue>Set</font>
#      <font color=blue>Description</font><br>
#    <p>
#      <font color=red>1</font> <A HREF="seasonal_mean_table.html">Tables</A> of DJF, MAM, JJA, SON and Annual Mean.<br>
#    <p>
#      <font color=red>2</font> <A HREF="annual_cycle.html">Line plots and Taylor diagrams</A> of Annual Cycle.<br>
#    <p>
#      <font color=red>3</font> <A HREF="annual_cycle_zt.html">Contour and Vertical profiles</A> of Annual Cycle.<br>
#    <p>
#      <font color=red>4</font> <A HREF="diurnal_cycle.html">Line and Harmonic Dail plots</A> of Diurnal Cycle.<br>
#    <p>
#      <font color=red>5</font> <A HREF="diurnal_cycle_zt.html">Contour plots</A> of Diurnal Cycle.<br>
#    <p>
#      <font color=red>6</font> <A HREF="pdf_daily.html">Line plots</A> of Probability Density Function.<br>
#    
#    </Table>
#
#    <p>
#    <p>
#    <Table>
#    <em>Click on Plot Type</em></b><p>
#      <A HREF="annual_cycle.html"><img src="../figures/pr_annual_cycle_sgp.png"  border=1 hspace=3 alt="Set 1" width="150" height="150"></a>
#      <A HREF="annual_cycle.html"><img src="../figures/pr_annual_cycle_taylor_diagram_sgp.png"  border=1 hspace=3 alt="Set 1" width="150" height="150"></a>
#      <A HREF="annual_cycle_zt.html"><img src="../figures/mod_cl_p_annual_cycle_clim_sgp.png"   border=1 hspace=3 alt="Set 3" width="150" height="150"></a>
#      <A HREF="DC_amip_line.html"><img src="../figures/ANN_cl_p_diff_sgp.png "   border=1 hspace=3 alt="Set 3" width="150" height="150"></a>
#      <A HREF="diurnal_cycle_zt.html"><img src="../figures/obs_cl_p_diurnal_clim_sgp.png"   border=1 hspace=3 alt="Set 3" width="150" height="150"></a>
#      <A HREF="pdf_daily.html"><img src="../figures/pr_JJA_pdf1_daily.png"   border=1 hspace=3 alt="Set 6" width="150" height="150"></a>
#    
#    </TH>
#    </TD>
#    
#    </TR>
#    </TABLE>
#    
#    </TD>
#    <Table>
#    <img src="../../misc/ARM_logo.png" >
#    </TABLE>
#    
#    
#    </html>"""
#
#
#
#    f.write(message)
#    f.close()
#    
#    return display_name, url
#
#
#
#
##    viewer = OutputViewer(path=root_dir)
##
##    # The name that's displayed on the viewer.
##    display_name = 'ARM Diags Annual Cycle'
##    set_name = 'arm_diags_annual_cycle'
##    cols = ['Description', 'Plot']
##    viewer.add_page(display_name, short_name=set_name, columns=cols)
##    viewer.add_group('Variable')
##
##    for param in parameter:
##        for var in param.variables:
##            viewer.add_row(var)
##            # Adding the description for this var to the current row.
##            # This was obtained and stored in the driver for this plotset.
##            viewer.add_col(param.viewer_descr[var])
##
##            ext = param.output_format[0]
##            file_name = os.path.join('..', set_name, param.case_id, '{}.{}'.format(var, ext))
##            viewer.add_col(file_name, is_file=True, title='Plot',
##                meta=create_metadata(param))
##
##    url = viewer.generate_page()
    
