import os
import collections
import cdms2
import cdutil
import acme_diags
from acme_diags.driver import utils
#from acme_diags.metrics import mean
import numpy as np
import json
from acme_diags.plot.cartopy import arm_diags_annual_cycle_plot

RefsTestMetrics = collections.namedtuple('RefsTestMetrics', ['refs', 'test', 'metrics'])

def rmse(predictions, targets):
    return np.sqrt(((predictions - targets) ** 2).mean())

def create_metrics(test, ref):
    """
    For this plotset, calculate the mean of the
    reference data and return a dict of that.
    """
    return {'test_mean': float(test.mean()),
            'ref_mean': float(ref.mean()),
            'test_std': float(test.std(ddof=1)),
            'ref_std': float(ref.std(ddof=1)),
            'rmse': float(rmse(test, ref)),
            'corr': float(np.corrcoef(test, ref)[0,1])
            }


def run_diag(parameter):
    variables = parameter.variables
    regions = parameter.regions
    ref_names = parameter.ref_names

    seasons = ['ANNUALCYCLE', 'ANN']
    # Both input data sets must be time-series files.
    # Raising an error will cause this specific set of
    # diagnostics with these parameters to be skipped.
    #if test_data.is_climo() or ref_data.is_climo():
    #    msg = 'Cannot run the plotset regional_mean_time_series '
    #    msg += 'because both the test and ref data need to be time-series files.'
    #    raise RuntimeError(msg)

    for region in regions:
        # The regions that are supported are in acme_diags/derivations/default_regions.py
        # You can add your own if it's not in there.
        print("Selected region: {}".format(region))
        vars_to_data = collections.OrderedDict()

        for season in seasons:
            print("Season: {}".format(season))
            for var in variables:
                print('Variable: {}'.format(var))
                test_data = utils.dataset.Dataset(parameter, test=True)
                #method 1, use built-in climo.py function to generate climatology, we know exactly about the algorithm
                test = test_data.get_climo_variable(var, season)
                print('test shape',test.shape, test.units)
                
                #method 2, use cdutil.ANNUALCYCLE.
                #test = test_data.get_timeseries_variable(var)
                ## Make sure data have correct montly Bounds
                #cdutil.setTimeBoundsMonthly(test)
                #print('test shape',test.shape, test.units)

                parameter.viewer_descr[var] = getattr(test, 'long_name', var)
                # Get the name of the data, appended with the years averaged.
                parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data)
                parameter.var_name = getattr(test, 'long_name', var)
                parameter.var_units = getattr(test, 'units', var)

                refs = []

                for ref_name in ref_names:    
                    setattr(parameter, 'ref_name', ref_name)
                    ref_data = utils.dataset.Dataset(parameter, ref=True)
                    print(ref_name)
                
                    parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data)
                    ref = ref_data.get_climo_variable(var, season)
                    ref.ref_name = ref_name

                    #ref = ref_data.get_timeseries_variable(var)
                    #cdutil.setTimeBoundsMonthly(ref)
                    #ref.ref_name = ref_name
                    
                    # TODO: Will this work if ref and test are timeseries data,
                    # but land_frac and ocean_frac are climo'ed.
                    print(ref_name)
                    test_domain = utils.general.select_point(region, test)
                    ref_domain = utils.general.select_point(region, ref)
                    refs.append(ref_domain)
                    #test_domain_year = cdutil.ANNUALCYCLE(test_domain)
                    #test_domain_year = test_domain
                    #ref_domain_year = cdutil.ANNUALCYCLE(ref_domain)
                    #ref_domain_year = ref_domain
                    #print(ref_domain)
                    #ref_domain_year.ref_name = ref_name

                    #refs.append(ref_domain_year)

                metrics_dict = create_metrics(test_domain,ref_domain)
                print(metrics_dict)
                #metrics_dict = ref_domain_year.mean()
                # print(test_domain_year.getTime().asComponentTime())
                # print(test.getTime().asComponentTime())

                result = RefsTestMetrics(test=test_domain, refs=refs, metrics=metrics_dict)
                vars_to_data[season] = result
                # Saving the metrics as a json.
                metrics_dict['unit'] = test.units
                parameter.output_file = '-'.join(
                            [ref_name, var, season, region])
                fnm = os.path.join(utils.general.get_output_dir(
                    parameter.current_set, parameter), parameter.output_file + '.json')
                with open(fnm, 'w') as outfile:
                     json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                # When running in a container, the paths are modified.
                fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
                    parameter, ignore_container=True), parameter.output_file + '.json')
                print('Metrics saved in: ' + fnm)
   
            #Compute and save stddev and correlation coefficient of models,for taylor diagram
            if season == 'ANNUALCYCLE':
                arm_diags_annual_cycle_plot.plot(var, vars_to_data[season], parameter)


        # TODO: How will this work when there are a bunch of plots for each image?
        # Yes, these files should be saved.
        # utils.general.save_ncfiles(parameter.current_set,
        #                     mv1_domain, mv2_domain, diff, parameter)
    return parameter


