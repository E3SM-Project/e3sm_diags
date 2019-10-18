from __future__ import print_function

import os
import math
import json
import numpy
import scipy.stats
import cdms2
import cdutil
import genutil
import acme_diags
from acme_diags.derivations import acme, default_regions
from acme_diags.driver import utils
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean, std
from acme_diags.plot.cartopy.enso_diags_plot import plot


def calculate_nino_index(test, nino_region_str, parameter):
    if test:
        data = utils.dataset.Dataset(parameter, test=True)
    else:
        data = utils.dataset.Dataset(parameter, ref=True)
    # type: cdms2.tvariable.TransientVariable
    sst = data.get_timeseries_variable('SST')
    nino_region = default_regions.regions_specs[nino_region_str]['domain']
    sst_nino = sst(nino_region)
    # Domain average
    sst_avg = cdutil.averager(sst_nino, axis='xy')
    # Get anomaly from annual cycle climatology
    sst_avg_anomaly = cdutil.ANNUALCYCLE.departures(sst_avg)
    nino_index = sst_avg_anomaly
    return data, sst, nino_index


def perform_regression(data, parameter, var, region, land_frac, ocean_frac, nino_index, test_sst, nino_region_str):
    ts_var = data.get_timeseries_variable(var)
    domain = utils.general.select_region(region, ts_var, land_frac, ocean_frac, parameter)
    # Average over selected region, and average
    # over months to get the yearly mean.
    cdutil.setTimeBoundsMonthly(domain)
    # Get anomaly from annual cycle climatology
    if parameter.print_statements:
        print('domain.shape: {}'.format(domain.shape))
    anomaly = cdutil.ANNUALCYCLE.departures(domain)
    nlat = len(anomaly.getLatitude())
    nlon = len(anomaly.getLongitude())
    reg_coe = anomaly[0, :, :](squeeze=1)
    confidence_levels = cdutil.ANNUALCYCLE.departures(domain)[0, :, :](squeeze=1)
    # Neither of the following methods work, so we just set values in confidence_levels
    # to be explicitly 0 or 1.
    #confidence_levels = anomaly[0, :, :](squeeze=1).fill(0) 
    #confidence_levels = numpy.zeros_like(reg_coe)
    for ilat in range(nlat):
        if parameter.print_statements:
            print('ilat: {}'.format(ilat))
        for ilon in range(nlon):
            dependent_var = anomaly[:, ilat, ilon]
            independent_var = nino_index
            # Uncomment the following line to use CDAT/genutil
            #slope, intercept = genutil.statistics.linearregression(dependent_var, x=independent_var)
            slope, _, _, pvalue, _ = scipy.stats.linregress(independent_var, dependent_var)
            reg_coe[ilat, ilon] = slope
            if parameter.print_statements:
                if slope > 5 and pvalue > 0.05:
                    print('slope, pvalue: {}, {}'.format(slope, pvalue))
            if pvalue < 0.05:
                # p-value < 5%
                # This implies significance at 95% confidence level
                #if parameter.print_statements:
                #    print('slope, pvalue: {}, {}'.format(slope, pvalue))
                confidence_levels[ilat, ilon] = 1
            else:
                confidence_levels[ilat, ilon] = 0
    if parameter.print_statements:
        print("confidence in fn:", confidence_levels.shape)
    # Use test.sst_units for both test and ref,
    # for consistent treatment of degC and C units
    reg_coe.units = '{}/{}'.format(ts_var.units, test_sst.units)
    if parameter.print_statements:
        print('reg_coe.shape: {}'.format(reg_coe.shape))
    # Create NetCDF files to save results
    file_name = '{}_{}_reg_coe.nc'.format(var, nino_region_str)
    with cdms2.open(file_name, 'w') as fout:
        fout.write(reg_coe)
    return domain, reg_coe, confidence_levels


def create_single_metrics_dict(values):
    d = {
        'min': float(min_cdms(values)),
        'max': float(max_cdms(values)),
        'mean': float(mean(values)),
        'std': float(std(values))
    }
    return d

def create_metrics(ref, test, ref_regrid, test_regrid, diff):
    """Creates the relevant metrics in a dictionary"""
    metrics_dict = {}
    metrics_dict['ref'] = create_single_metrics_dict(ref)
    metrics_dict['ref_regrid'] = create_single_metrics_dict(ref_regrid)
    metrics_dict['test'] = create_single_metrics_dict(test)
    metrics_dict['test_regrid'] = create_single_metrics_dict(test_regrid)
    metrics_dict['diff'] = create_single_metrics_dict(diff)
    d = metrics_dict['diff']
    d['rmse'] = float(rmse(test_regrid, ref_regrid))
    mean_sq = d['mean'] ** 2
    rmse_sq = d['rmse'] ** 2
    d['fraction'] = mean_sq / rmse_sq
    return metrics_dict


def run_diag(parameter):
    variables = parameter.variables
    seasons = parameter.seasons
    regions = parameter.regions
    nino_region_str = parameter.nino_region

    # Test
    test_data, test_sst, test_nino_index = calculate_nino_index(
        True, nino_region_str, parameter)

    # Reference
    ref_data, _, ref_nino_index = calculate_nino_index(
        False, nino_region_str, parameter)

    for season in seasons:
        if parameter.print_statements:
            print('Season: {}'.format(season))
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

        for var in variables:
            if parameter.print_statements:
                print('Variable: {}'.format(var))
            parameter.var_id = var

            parameter.viewer_descr[var] = 'Regression coefficient, {} over nino index'.format(var)

            for region in regions:
                if parameter.print_statements:
                    print("Selected region: {}".format(region))

                # This will be the title of the plot.
                parameter.main_title = ' '.join(
                    [parameter.viewer_descr[var], region])

                # Test
                test_domain, test_reg_coe, test_confidence_levels  = perform_regression(test_data, parameter, var, region, land_frac, ocean_frac, test_nino_index, test_sst, nino_region_str)

                # Reference
                ref_domain, ref_reg_coe, ref_confidence_levels = perform_regression(ref_data, parameter, var, region, land_frac, ocean_frac, ref_nino_index, test_sst, nino_region_str)

                # Difference
                # Regrid towards the lower resolution of the two variables for calculating the difference.
                test_reg_coe_regrid, ref_reg_coe_regrid = utils.general.regrid_to_lower_res(test_reg_coe, ref_reg_coe, parameter.regrid_tool, parameter.regrid_method)
                diff = test_reg_coe_regrid - ref_reg_coe_regrid

                # Metrics
                metrics_dict = create_metrics(ref_reg_coe, test_reg_coe, ref_reg_coe_regrid, test_reg_coe_regrid, diff)
                # If nor defined, determined contour_levels
                if not parameter.contour_levels:
                    # We want contour levels for the plot,
                    # which uses original (non-regridded) test and ref,
                    # so we use those min and max values.
                    min_contour_level = math.floor(min(
                        metrics_dict['ref']['min'],
                        metrics_dict['test']['min']
                    ))
                    max_contour_level = math.ceil(max(
                        metrics_dict['ref']['max'],
                        metrics_dict['test']['max']
		    ))
                    contour_level_range = max_contour_level - min_contour_level
                    step_size = (contour_level_range // 15) + 1
                    parameter.contour_levels = list(range(min_contour_level, max_contour_level + 1, step_size))
                parameter.output_file = 'regression-coefficient-{}-over-nino-index'.format(var.lower())
                
                # Saving the metrics as a json.
                metrics_dict['unit'] = test_reg_coe_regrid.units
                metrics_output_file_name = os.path.join(utils.general.get_output_dir(parameter.current_set, parameter), parameter.output_file + '.json')
                with open(metrics_output_file_name, 'w') as outfile:
                    json.dump(metrics_dict, outfile)
                # Get the filename that the user has passed in and display that.
                # When running in a container, the paths are modified.
                metrics_output_file_name = os.path.join(utils.general.get_output_dir(parameter.current_set, parameter, ignore_container=True), parameter.output_file + '.json')
                print('Metrics saved in: {}'.format(metrics_output_file_name))

                # Plot
                parameter.var_region = region
                # Plot original ref and test, not regridded versions.
                plot(ref_reg_coe, test_reg_coe, diff,
                     metrics_dict, ref_confidence_levels, test_confidence_levels,
                     parameter)
                utils.general.save_ncfiles(parameter.current_set,
                                           test_reg_coe, ref_reg_coe, diff, parameter)
    return parameter
