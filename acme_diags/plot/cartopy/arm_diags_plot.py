import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp2d
from acme_diags.driver.utils.general import get_output_dir
from acme_diags.plot.cartopy.taylor_diagram import TaylorDiagram
from acme_diags.driver.utils.diurnal_cycle import fastAllGridFT



def plot_annual_cycle(var, vars_to_data, parameter):
    line_color = ['r', 'b', 'g', 'm']

    num_year = int(parameter.test_end_yr) - int(parameter.test_start_yr) +1
    fig = plt.figure()# Create figure
    ax1  =fig.add_axes([0.15, 0.14, 0.8, 0.8]) # Create axes
    xax =  np.arange (1,13,1)

    print(vars_to_data)
    refs = vars_to_data.refs
    test = vars_to_data.test
    ax1.plot(xax, test.asma(), 'k', linewidth=2,label = 'model' +' ({0:.1f})'.format(np.mean(test.asma())))
    for i_ref, ref in enumerate(refs):
        ax1.plot(xax, ref.asma(), line_color[i_ref], linewidth=2,label = ref.ref_name +' ({0:.1f})'.format(np.mean(ref.asma())))
    my_xticks = ['J','F','M','A','M','J','J','A','S','O','N','D']
    plt.xticks(xax, my_xticks)
    plt.xlim(1,12)
#     plt.ylim(ylim[va_ind])
    plt.title('Annual Cycle: Model vs OBS' )
    plt.xlabel('Month')
    plt.legend(loc='best',prop={'size':10})
    plt.ylabel(parameter.var_name + ' (' +parameter.var_units+ ')')


    # Save the figure.
    output_file_name = parameter.output_file
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set,
            parameter), output_file_name + '.' + f)
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
            ignore_container=True), output_file_name + '.' + f)
        print('Plot saved in: ' + fnm)

    plt.close()

def plot_diurnal_cycle(var, vars_to_data, parameter):
    test = vars_to_data.test[0]
    ref = vars_to_data.refs[0][0]
    lst = vars_to_data.misc[0]
    t_conv = lst[0][0]
    print(ref)
    print(test)
    print(lst)
    print(t_conv)
    c, maxvalue, tmax = fastAllGridFT(test,lst)
    print(c,maxvalue, tmax)
    c, maxvalue, tmax = fastAllGridFT(ref,lst)
    print(c,maxvalue, tmax)
    output_file_name = parameter.output_file+'-'+'diurnal-cycle'

    fig = plt.figure()# Create figure
    ax  =fig.add_axes([0.15, 0.14, 0.8, 0.8]) # Create axes
 
    for index in range(2):
       if index == 0:
           data = test
           line_c = 'k'
           data_name = parameter.test_name
       else:
           data = ref
           line_c = 'r'
           data_name = parameter.ref_name
          
       time_freq = len(data)
       c, maxvalue, tmax = fastAllGridFT(data,[0])
       xax = np.linspace(0,48,time_freq *2)
       ax.plot(xax,np.concatenate((data,data)),'.'+line_c, label = data_name)
       xax = np.linspace(0,48,time_freq *2*3)
       print('tmax',tmax[0])
       w = 2.0*np.pi/24
       yax = (c + maxvalue[0] *np.sin(w*xax+np.pi/2-tmax[0]*w))[0]
       print(yax)
       ax.plot(xax,yax, line_c, label = 'First harmonic')
       plt.xlim([24-t_conv,47-t_conv+1])
       plt.ylim([-0.5,7])
       plt.xlabel('local solar time [hr]')
       #plt.ylabel(parameter.var_name + ' (' +parameter.var_units+ ')')
       plt.ylabel('Total Precipitation Rate' + ' (' +parameter.var_units+ ')')
       xax = np.arange(24-t_conv,47-t_conv,3)
       my_xticks = ['0h','3h','6h','9h','12h','15h','18h','21h']
       plt.xticks(xax, my_xticks)
       plt.legend(loc="upper right")
       plt.title(output_file_name.replace('-', ' '))
    

    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        fnm = os.path.join(get_output_dir(parameter.current_set,
            parameter), output_file_name + '.' + f)
        plt.savefig(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
            ignore_container=True), output_file_name + '.' + f)
        print('Plot saved in: ' + fnm)

    plt.close()

    

    
   
    


def plot_diurnal_cycle_zt(var, vars_to_data, parameter):
    ref = vars_to_data.refs[0]
    test = vars_to_data.test
    lst = vars_to_data.misc
    month = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

    for index in range(2):
        fig, axs = plt.subplots(4,3, figsize=(15, 12), facecolor='w', edgecolor='k',sharex=True,sharey=True)
        fig.subplots_adjust(hspace = .3, wspace=.1)
        axs = axs.ravel()
        t_conv = lst[0][0][0]
        print(t_conv)
        for imon in range(12):
            if index==0:
                 title= parameter.ref_name
                 data = ref
                 data_name = 'ref'
                 if 'armdiags' in title: 
                     data = data[:,:,::-1]

            else:
                 title= parameter.test_name
                 data = test
                 data_name = 'test'

            time_freq = data.shape[1]
            yy=np.linspace(0,48,time_freq *2)
            xx=np.linspace(100,1000,37)
            x,y=np.meshgrid(xx,yy)
            data_con=np.concatenate((data[imon,:,:],data[imon,:,:]),axis=0)
            im=axs[imon].pcolormesh(y,x,data_con[:,:], vmin=0, vmax=25,cmap='jet', shading='auto')
            axs[imon].set_title(month[imon])
            plt.xlim([24-t_conv,47-t_conv]) 
            xax = np.arange(24-t_conv,47-t_conv,3)
            my_xticks = ['0','3','6','9','12','15','18','21']
            plt.xticks(xax, my_xticks)
            #plt.setp(axs[imon].get_xticklabels(), visible=True)

        for ax in axs[9:12]:
            ax.set_xlabel('Local time (hr)')
        for ax in axs[::3]:
            ax.set_ylabel('Pressure (mb)')
        axs[0].invert_yaxis()
        plt.suptitle(title)
        fig.subplots_adjust(right=0.8)
        cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
        fig.colorbar(im, cax=cbar_ax)
        plt.title('cl (%)')

        output_file_name = parameter.output_file+'-'+data_name
        for f in parameter.output_format:
            f = f.lower().split('.')[-1]
            fnm = os.path.join(get_output_dir(parameter.current_set,
                parameter), output_file_name + '.' + f)
            plt.savefig(fnm)
            # Get the filename that the user has passed in and display that.
            # When running in a container, the paths are modified.
            fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
                ignore_container=True), output_file_name + '.' + f)
            print('Plot saved in: ' + fnm)
    
        plt.close()

