from __future__ import annotations

import csv
import os
import xarray as xr
import numpy
import matplotlib.pyplot as plt
from scipy.stats import binned_statistic
import scipy.io
import math
import json
from typing import TYPE_CHECKING  # , Optional

from e3sm_diags.driver import utils

if TYPE_CHECKING:
    from e3sm_diags.parameter.mp_partition_parameter import MPpartitionParameter

import cdutil
import numpy

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)

# This analysis set for mixed-phase cloud partition/T5050 metrics is requested by the E3SM Aerosol Working Group. The script is integrated in e3sm_diags by Jill Zhang and Yuying Zhang, with contribution from Yunpeng Shan, Jiwen Fan, Xue Zheng and Susannah Burrows.

def flatten_array(var):
    var_1d = var.stack(stacked=[...]).values
    var_1d = var_1d[~numpy.isnan(var_1d)]
    return var_1d


def compute_lcf(cice, cliq, temp, landfrac):
    ctot = cice + cliq
    ctot_sel = ctot.where((temp >= 220) & (temp <= 280)).where(ctot > 1e-9).where(landfrac==0)
    cliq_sel = cliq.where(ctot_sel.notnull())
    temp_sel = temp.where(ctot_sel.notnull())
    
    ctot_1d = flatten_array(ctot_sel)
    cliq_1d = flatten_array(cliq_sel)
    temp_1d = flatten_array(temp_sel)
    
    lcf = cliq_1d/ctot_1d
    
    mean_stat = binned_statistic(temp_1d, lcf, 
                             statistic='mean', 
                             bins=20, 
                             range=(220, 280))
    
    temp_bin_center = (mean_stat.bin_edges[:-1] + mean_stat.bin_edges[1:]) / 2
    
    return temp_bin_center, mean_stat.statistic


def run_diag(parameter: MPpartitionParameter) -> MPpartitionParameter:
    """Runs the mixed-phase partition/T5050 diagnostic.

    :param parameter: Parameters for the run
    :type parameter: CoreParameter
    :raises ValueError: Invalid run type
    :return: Parameters for the run
    :rtype: CoreParameter
    """
    test_data_path = parameter.test_data_path
    reference_data_path = parameter.reference_data_path
    #variables = parameter.variables[0].split(", ")
    run_type = parameter.run_type

    metrics_dict_test = {}
    metrics_dict_ref = {}
    test_data = utils.dataset.Dataset(parameter, test=True)
    test = test_data.get_timeseries_variable("T")
    parameter.test_name_yrs = utils.general.get_name_and_yrs(
        parameter, test_data, season
    )
    parameter.ref_name_yrs = "OBS"

    datapath = '/Users/zhang40/Documents/e3sm_diags_data/e3sm_diags_test_data/E3SM_v2/'
    landfrac = xr.open_dataset(f'{datapath}/LANDFRAC_200001_200012.nc').sel(lat=slice(-70,-30))['LANDFRAC']
    temp = xr.open_dataset(f'{datapath}/T_2000001_200012.nc').sel(lat=slice(-70,-30))['T']
    cice = xr.open_dataset(f'{datapath}/CLDICE_2000001_2000012.nc').sel(lat=slice(-70,-30))['CLDICE']
    cliq = xr.open_dataset(f'{datapath}/CLDLIQ_2000001_2000012.nc').sel(lat=slice(-70,-30))['CLDLIQ']

    metrics_dir = 
#
#    for aerosol in variables:
#        logger.info("Variable: {}".format(aerosol))
#        metrics_dict_test[aerosol] = generate_metrics_dic(
#            test_data, aerosol, season
#        )
#
#    if run_type == "model_vs_model":
#        ref_data = utils.dataset.Dataset(parameter, ref=True)
#        parameter.ref_name_yrs = utils.general.get_name_and_yrs(
#            parameter, ref_data, season
#        )
#        for aerosol in variables:
#            metrics_dict_ref[aerosol] = generate_metrics_dic(
#                ref_data, aerosol, season
#            )
#
#    elif run_type == "model_vs_obs":
#        for aerosol in variables:
#            metrics_dict_ref[aerosol] = {
#                "Surface Emission (Tg/yr)": f"{MISSING_VALUE:.3f}",
#                "Elevated Emission (Tg/yr)": f"{MISSING_VALUE:.3f}",
#                "Sink (Tg/s)": f"{MISSING_VALUE:.3f}",
#                "Dry Deposition (Tg/yr)": f"{MISSING_VALUE:.3f}",
#                "Wet Deposition (Tg/yr)": f"{MISSING_VALUE:.3f}",
#                "Burden (Tg)": f"{MISSING_VALUE:.3f}",
#                "Lifetime (Days)": f"{MISSING_VALUE:.3f}",
#            }
#    else:
#        raise ValueError("Invalid run_type={}".format(run_type))
#
#    parameter.output_file = f"{parameter.test_name}-{season}-budget-table"
#    fnm = os.path.join(
#        utils.general.get_output_dir(parameter.current_set, parameter),
#        parameter.output_file + ".csv",
#    )
#
#    with open(fnm, "w") as table_csv:
#        writer = csv.writer(
#            table_csv,
#            delimiter=",",
#            quotechar="'",
#            quoting=csv.QUOTE_MINIMAL,
#            lineterminator="\n",
#        )
#        writer.writerow([f"Test: {parameter.test_name_yrs}"])
#        writer.writerow([f"Ref: {parameter.ref_name_yrs}"])
#        writer.writerow(
#            [
#                " ",
#                "test",
#                "ref",
#            ]
#        )
#        for key, values in metrics_dict_test.items():
#            writer.writerow([SPECIES_NAMES[key]])
#            for value in values:
#                line = []
#                line.append(value)
#                line.append(values[value])
#                line.append(metrics_dict_ref[key][value])
#                writer.writerows([line])
#            writer.writerows([""])
#
#    logger.info(f"Metrics saved in {fnm}")

    return parameter
