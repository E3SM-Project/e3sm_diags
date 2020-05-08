from __future__ import print_function

import os
import json
import cdms2
import MV2
import acme_diags
from acme_diags.plot import plot
from acme_diags.derivations import acme
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean, std
from acme_diags.driver import utils
import collections

Results_Container = collections.namedtuple('Results_Container', ['refs', 'tests', 'diff', 'metrics'])


def create_metrics(ref_regrid, test_regrid, diff):
    """Creates the mean, max, min, rmse, corr in a dictionary"""
    metrics_dict = {}
    metrics_dict['ref_regrid'] = {
        'min': float(min_cdms(ref_regrid)),
        'max': float(max_cdms(ref_regrid)),
        'mean': float(mean(ref_regrid)),
        'std': float(std(ref_regrid))
    }
    metrics_dict['test_regrid'] = {
        'min': float(min_cdms(test_regrid)),
        'max': float(max_cdms(test_regrid)),
        'mean': float(mean(test_regrid)),
        'std': float(std(test_regrid))
    }
    metrics_dict['diff'] = {
        'min': float(min_cdms(diff)),
        'max': float(max_cdms(diff)),
        'mean': float(mean(diff))
    }
    metrics_dict['misc'] = {
        'rmse': float(rmse(test_regrid, ref_regrid)),
        'corr': float(corr(test_regrid, ref_regrid))
    }
    return metrics_dict


def run_diag(parameter):
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, 'ref_name', '')
    regions = parameter.regions

    test_data = utils.dataset.Dataset(parameter, test=True)
    ref_data = utils.dataset.Dataset(parameter, ref=True)

    for season in seasons:
        # Get the name of the data, appended with the years averaged.
        parameter.test_name_yrs = utils.general.get_name_and_yrs(parameter, test_data, season)
        parameter.ref_name_yrs = utils.general.get_name_and_yrs(parameter, ref_data, season)

        # Get land/ocean fraction for masking.
        try:
            land_frac = test_data.get_climo_variable('LANDFRAC', season)
            ocean_frac = test_data.get_climo_variable('OCNFRAC', season)
        except:
            mask_path = os.path.join(acme_diags.INSTALL_PATH, 'acme_ne30_ocean_land_mask.nc')
            with cdms2.open(mask_path) as f:
                land_frac = f('LANDFRAC')
                ocean_frac = f('OCNFRAC')

        tests = []
        refs = []
        ZONAL = 0
        MERIDIONAL = 1
        var_title = '{}-{}'.format(parameter.variables[ZONAL], parameter.variables[MERIDIONAL])
        for var in variables:
            print('Variable: {}'.format(var))
            parameter.var_id = var

            mv1 = test_data.get_climo_variable(var, season)
            mv2 = ref_data.get_climo_variable(var, season)

            parameter.viewer_descr[var] = mv1.long_name if hasattr(
                mv1, 'long_name') else 'No long_name attr in test data.'
            s = parameter.viewer_descr[var].replace('Zonal ', '').replace('Meridional ', '').capitalize()
            parameter.viewer_descr[var] = s


            # For variables with a z-axis.
            if mv1.getLevel() and mv2.getLevel():
                plev = parameter.plevs
                print('Selected pressure level: {}'.format(plev))

                mv1_p = utils.general.convert_to_pressure_levels(mv1, plev, test_data, var, season)
                mv2_p = utils.general.convert_to_pressure_levels(mv2, plev, ref_data, var, season)

                # Select plev.
                for ilev in range(len(plev)):
                    mv1 = mv1_p[ilev, ]
                    mv2 = mv2_p[ilev, ]

                    for region in regions:
                        print("Selected region: {}".format(region))

                        mv1_domain = utils.general.select_region(region, mv1, land_frac, ocean_frac, parameter)
                        mv2_domain = utils.general.select_region(region, mv2, land_frac, ocean_frac, parameter)

                        parameter.output_file = '-'.join(
                            [ref_name, var_title, str(int(plev[ilev])), season, region])
                        #parameter.main_title = str(
                        #    ' '.join([var, str(int(plev[ilev])), 'mb', season, region]))

                        if var.find('TAU') == 0:
                            parameter.main_title = str(' '.join(['TAU', str(int(plev[ilev])), 'mb', season, region]))
                        elif var.find('U') ==0 or var.find('V') ==0:
                            parameter.main_title = str(' '.join(['Wind', str(int(plev[ilev])), 'mb', season, region]))
                        else:
                            parameter.main_title = str(' '.join([var, str(int(plev[ilev])), 'mb', season, region]))
    
                        # Regrid towards the lower resolution of the two
                        # variables for calculating the difference.
                        mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                            mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                        tests.append(mv1_reg)
                        refs.append(mv2_reg)


            # For variables without a z-axis.
            else:
                for region in regions:
                    print("Selected region: {}".format(region))

                    mv1_domain = utils.general.select_region(region, mv1, land_frac, ocean_frac, parameter)
                    mv2_domain = utils.general.select_region(region, mv2, land_frac, ocean_frac, parameter)

                    parameter.output_file = '-'.join(
                        [ref_name, var_title, season, region])
                    #parameter.main_title = str(' '.join([var, season, region]))
                    if var.find('TAU') == 0:
                        parameter.main_title = str(' '.join(['TAU', season, region]))
                    elif var.find('U') == 0 or var.find('V') == 0:
                        parameter.main_title = str(' '.join(['Near surface wind', season, region]))
                    else:
                        parameter.main_title = str(' '.join([var, season, region]))
  

                    # Regrid towards the lower resolution of the two
                    # variables for calculating the difference.
                    mv1_reg, mv2_reg = utils.general.regrid_to_lower_res(
                        mv1_domain, mv2_domain, parameter.regrid_tool, parameter.regrid_method)

                    tests.append(mv1_reg)
                    refs.append(mv2_reg)

        # Euclidean norm of the vector
        test = (tests[ZONAL] **2.0 + tests[MERIDIONAL]**2.0)**0.5
        ref = (refs[ZONAL] **2.0 + refs[MERIDIONAL]**2.0)**0.5
        diffs = []
        diff = test - ref
        metrics_dict = create_metrics(ref, test, diff)
        tests.append(test)
        refs.append(ref)
        diffs.append(tests[ZONAL] - ref[ZONAL])
        diffs.append(tests[MERIDIONAL] - ref[MERIDIONAL])
        # Note that diff is not the vector norm of the above two components.
        # It is the difference of the test and ref vector norms.
        diffs.append(diff)
        #result = Results_Container(tests=tests, refs=refs, diff=diff, metrics = metrics_dict)

        # Saving the metrics as a json.
        metrics_dict['unit'] = tests[0].units
        fnm = os.path.join(utils.general.get_output_dir(
            parameter.current_set, parameter), parameter.output_file + '.json')
        with open(fnm, 'w') as outfile:
                json.dump(metrics_dict, outfile)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(utils.general.get_output_dir(parameter.current_set,
            parameter, ignore_container=True), parameter.output_file + '.json')
        print('Metrics saved in: ' + fnm)

        parameter.var_region = region
        # Calling plot:
        # Call acme_diags.__init__/plot, which uses
        # acme_diags.__init__/_get_plot_fcn to call
        # acme_diags.plot.{backend}.{set_name}_plot
        # {parameter.backend} defaults to `mpl`
        # (in acme_diags.parameter.core_parameter.CoreParameter.__init__),
        # which will use the `cartopy` backend.
        # `set_name` is set to `parameter.curent_set` below.
        # `parameter.current_set` is set to the current set from
        # `parameters.sets` in acme_diags.acme_diags_driver.run_diag
        # e.g. acme_diags.plot.cartopy.lat_lon_vector_plot is called
        plot(parameter.current_set, tests,
             refs, diffs, metrics_dict, parameter)
        utils.general.save_ncfiles(parameter.current_set,
                          test, ref, diff, parameter)
     
    return parameter
