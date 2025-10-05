"""
This analysis set for mixed-phase cloud partition/T5050 metrics is requested by
the E3SM Aerosol Working Group. The script is integrated in e3sm_diags by Jill
Zhang and Yuying Zhang, with contribution from Yunpeng Shan, Jiwen Fan,
Xue Zheng and Susannah Burrows.
"""

from __future__ import annotations

import glob
import json
import os
from typing import TYPE_CHECKING

import numpy
import xarray as xr
from scipy.stats import binned_statistic

from e3sm_diags import INSTALL_PATH
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.mp_partition_plot import plot

if TYPE_CHECKING:
    from e3sm_diags.parameter.mp_partition_parameter import MPpartitionParameter


logger = _setup_child_logger(__name__)


def flatten_array(var):
    var_1d = var.stack(stacked=[...]).values
    var_1d = var_1d[~numpy.isnan(var_1d)]
    return var_1d


def compute_lcf(cice, cliq, temp, landfrac):
    ctot = cice + cliq
    ctot_sel = (
        ctot.where((temp >= 220) & (temp <= 280))
        .where(ctot > 1e-9)
        .where(landfrac == 0)
    )
    cliq_sel = cliq.where(ctot_sel.notnull())
    temp_sel = temp.where(ctot_sel.notnull())

    ctot_1d = flatten_array(ctot_sel)
    cliq_1d = flatten_array(cliq_sel)
    temp_1d = flatten_array(temp_sel)

    lcf = cliq_1d / ctot_1d

    mean_stat = binned_statistic(
        temp_1d, lcf, statistic="mean", bins=20, range=(220, 280)
    )

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
    run_type = parameter.run_type
    season = "ANN"

    # Read reference data first

    benchmark_data_path = os.path.join(
        INSTALL_PATH,
        "control_runs",
        "mixed-phase_partition_data_1985-2014.json",
    )

    with open(benchmark_data_path, "r") as myfile:
        lcf_file = myfile.read()

    # parse file
    metrics_dict = json.loads(lcf_file)

    test_data = Dataset(parameter, data_type="test")
    # test = test_data.get_timeseries_variable("LANDFRAC")
    # print(dir(test))
    # landfrac = test_data.get_timeseries_variable("LANDFRAC")(cdutil.region.domain(latitude=(-70.0, -30, "ccb")))
    # temp = test_data.get_timeseries_variable("T")(cdutil.region.domain(latitude=(-70.0, -30, "ccb")))
    # cice = test_data.get_timeseries_variable("CLDICE")(cdutil.region.domain(latitude=(-70.0, -30, "ccb")))
    # cliq = test_data.get_timeseries_variable("CLDLIQ")(cdutil.region.domain(latitude=(-70.0, -30, "ccb")))

    test_data_path = parameter.test_data_path
    start_year = parameter.test_start_yr
    end_year = parameter.test_end_yr
    # TODO the time subsetting and variable derivation should be replaced during cdat revamp
    try:
        # xr.open_mfdataset() can accept an explicit list of files.
        landfrac = xr.open_mfdataset(glob.glob(f"{test_data_path}/LANDFRAC_*")).sel(
            lat=slice(-70, -30), time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
        )["LANDFRAC"]
        temp = xr.open_mfdataset(glob.glob(f"{test_data_path}/T_*.nc")).sel(
            lat=slice(-70, -30), time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
        )["T"]
        cice = xr.open_mfdataset(glob.glob(f"{test_data_path}/CLDICE_*.nc")).sel(
            lat=slice(-70, -30), time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
        )["CLDICE"]
        cliq = xr.open_mfdataset(glob.glob(f"{test_data_path}/CLDLIQ_*.nc")).sel(
            lat=slice(-70, -30), time=slice(f"{start_year}-01-01", f"{end_year}-12-31")
        )["CLDLIQ"]
    except OSError:
        logger.info(
            f"No files to open for variables within {start_year} and {end_year} from {test_data_path}."
        )
        raise

    parameter.test_name_yrs = test_data.get_name_yrs_attr(season)

    metrics_dict["test"] = {}
    metrics_dict["test"]["T"], metrics_dict["test"]["LCF"] = compute_lcf(
        cice, cliq, temp, landfrac
    )

    if run_type == "model-vs-model":
        ref_data = Dataset(parameter, data_type="ref")

        ref_data_path = parameter.reference_data_path
        start_year = parameter.ref_start_yr
        end_year = parameter.ref_end_yr
        # xr.open_mfdataset() can accept an explicit list of files.
        try:
            landfrac = xr.open_mfdataset(glob.glob(f"{ref_data_path}/LANDFRAC_*")).sel(
                lat=slice(-70, -30),
                time=slice(f"{start_year}-01-01", f"{end_year}-12-31"),
            )["LANDFRAC"]
            temp = xr.open_mfdataset(glob.glob(f"{ref_data_path}/T_*.nc")).sel(
                lat=slice(-70, -30),
                time=slice(f"{start_year}-01-01", f"{end_year}-12-31"),
            )["T"]
            cice = xr.open_mfdataset(glob.glob(f"{ref_data_path}/CLDICE_*.nc")).sel(
                lat=slice(-70, -30),
                time=slice(f"{start_year}-01-01", f"{end_year}-12-31"),
            )["CLDICE"]
            cliq = xr.open_mfdataset(glob.glob(f"{ref_data_path}/CLDLIQ_*.nc")).sel(
                lat=slice(-70, -30),
                time=slice(f"{start_year}-01-01", f"{end_year}-12-31"),
            )["CLDLIQ"]
        except OSError:
            logger.info(
                f"No files to open for variables within {start_year} and {end_year} from {ref_data_path}."
            )
            raise

        # landfrac = ref_data.get_timeseries_variable("LANDFRAC")(
        #    cdutil.region.domain(latitude=(-70.0, -30, "ccb"))
        # )
        # temp = ref_data.get_timeseries_variable("T")(
        #    cdutil.region.domain(latitude=(-70.0, -30, "ccb"))
        # )
        # cice = ref_data.get_timeseries_variable("CLDICE")(
        #    cdutil.region.domain(latitude=(-70.0, -30, "ccb"))
        # )
        # cliq = ref_data.get_timeseries_variable("CLDLIQ")(
        #    cdutil.region.domain(latitude=(-70.0, -30, "ccb"))
        # )
        parameter.ref_name_yrs = ref_data.get_name_yrs_attr(season)
        metrics_dict["ref"] = {}
        metrics_dict["ref"]["T"], metrics_dict["ref"]["LCF"] = compute_lcf(
            cice, cliq, temp, landfrac
        )
    parameter.output_file = "mixed-phase_partition"

    # TODO: save metrics
    plot(metrics_dict, parameter)

    return parameter
