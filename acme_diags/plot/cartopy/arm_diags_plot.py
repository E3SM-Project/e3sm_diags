import os
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import numpy as np
from acme_diags.driver.utils.general import get_output_dir


def plot(var, vars_to_data, parameter):
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


#        x = np.arange(num_year)
#        ax1.set_xticks(x)
#        x_ticks_labels = [str(x) for x in range(int(parameter.test_start_yr),int(parameter.test_end_yr)+1)]
#        ax1.set_xticklabels(x_ticks_labels, rotation='45', fontsize=6)
#
#    # Figure title.
#    fig.suptitle('Annual mean ' + var + ' over regions ' + parameter.test_name_yrs, x=0.5, y=0.97, fontsize=15)

    # Save the figure.
    output_file_name = var
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

