"""
This analysis set for mixed-phase cloud partition/T5050 metrics is requested by
the E3SM Aerosol Working Group. The script is integrated in e3sm_diags by Jill
Zhang and Yuying Zhang, with contribution from Yunpeng Shan, Jiwen Fan,
Xue Zheng and Susannah Burrows.
"""

from __future__ import annotations

import json
import os
from typing import TYPE_CHECKING

import numpy
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
    """
    Compute Liquid Condensate Fraction (LCF).

    TODO: Implement IDL algorithm changes when using COSP variables:
    - Use threshold 0.001 instead of 1e-9 for COSP variables
    - Implement manual binning following IDL algorithm:
      tind = int((temp - 220) / 3)
      Accumulate lcf values per bin, then divide by count
    """
    ctot = cice + cliq
    ctot_sel = (
        ctot.where((temp >= 220) & (temp <= 280))
        .where(ctot > 1e-9)  # TODO: Use 0.001 for COSP variables
        .where(landfrac == 0)
    )
    cliq_sel = cliq.where(ctot_sel.notnull())
    temp_sel = temp.where(ctot_sel.notnull())

    ctot_1d = flatten_array(ctot_sel)
    cliq_1d = flatten_array(cliq_sel)
    temp_1d = flatten_array(temp_sel)

    lcf = cliq_1d / ctot_1d

    # TODO: Replace with manual binning for COSP variables following IDL algorithm
    mean_stat = binned_statistic(
        temp_1d, lcf, statistic="mean", bins=20, range=(220, 280)
    )

    temp_bin_center = (mean_stat.bin_edges[:-1] + mean_stat.bin_edges[1:]) / 2

    return temp_bin_center, mean_stat.statistic


def run_diag(parameter: MPpartitionParameter) -> MPpartitionParameter:
    """
    Runs the mixed-phase partition/T5050 diagnostic.

    Parameters
    ----------
    parameter : CoreParameter
        Parameters for the run.

    Raises
    ------
    ValueError
        If the run type is invalid.

    Returns
    -------
    CoreParameter
        Parameters for the run.
    """
    run_type = parameter.run_type
    season = "ANN"

    benchmark_data_path = os.path.join(
        INSTALL_PATH,
        "control_runs",
        "mixed-phase_partition_data_1985-2014.json",
    )

    with open(benchmark_data_path, "r") as myfile:
        lcf_file = myfile.read()

    metrics_dict = json.loads(lcf_file)

    test_data = Dataset(parameter, data_type="test")

    # Load required variables using Dataset class
    # Try new COSP variables first, fall back to original variables
    use_cosp_variables = False

    try:
        landfrac_ds = test_data.get_time_series_dataset("LANDFRAC", single_point=True)
        landfrac = landfrac_ds["LANDFRAC"].sel(lat=slice(-70, -30))
    except Exception as e:
        logger.warning(f"Could not load LANDFRAC: {e}")
        landfrac = None

    try:
        # Try CLD_CAL_TMPICE first, fall back to CLDICE
        try:
            cice_ds = test_data.get_time_series_dataset(
                "CLD_CAL_TMPICE", single_point=True
            )
            cice = cice_ds["CLD_CAL_TMPICE"].sel(lat=slice(-70, -30))
            use_cosp_variables = True
            logger.info("Using CLD_CAL_TMPICE for cloud ice")
        except Exception:
            cice_ds = test_data.get_time_series_dataset("CLDICE", single_point=True)
            cice = cice_ds["CLDICE"].sel(lat=slice(-70, -30))
            logger.info("Using CLDICE for cloud ice")

        # Try CLD_CAL_TMPLIQ first, fall back to CLDLIQ
        try:
            cliq_ds = test_data.get_time_series_dataset(
                "CLD_CAL_TMPLIQ", single_point=True
            )
            cliq = cliq_ds["CLD_CAL_TMPLIQ"].sel(lat=slice(-70, -30))
            use_cosp_variables = True
            logger.info("Using CLD_CAL_TMPLIQ for cloud liquid")
        except Exception:
            cliq_ds = test_data.get_time_series_dataset("CLDLIQ", single_point=True)
            cliq = cliq_ds["CLDLIQ"].sel(lat=slice(-70, -30))
            logger.info("Using CLDLIQ for cloud liquid")

        # Load temperature: cosp_temp for COSP variables, T for original variables
        if use_cosp_variables:
            try:
                temp_ds = test_data.get_time_series_dataset(
                    "cosp_temp", single_point=True
                )
                temp = temp_ds["cosp_temp"].sel(lat=slice(-70, -30))
                logger.info("Using cosp_temp for temperature")
            except Exception:
                logger.warning("Could not load cosp_temp, falling back to T")
                temp_ds = test_data.get_time_series_dataset("T", single_point=True)
                temp = temp_ds["T"].sel(lat=slice(-70, -30))
        else:
            temp_ds = test_data.get_time_series_dataset("T", single_point=True)
            temp = temp_ds["T"].sel(lat=slice(-70, -30))
            logger.info("Using T for temperature")

    except Exception as e:
        logger.info(f"Error loading variables for test data: {e}")
        raise

    parameter.test_name_yrs = test_data.get_name_yrs_attr(season)

    metrics_dict["test"] = {}
    metrics_dict["test"]["T"], metrics_dict["test"]["LCF"] = compute_lcf(
        cice, cliq, temp, landfrac
    )

    if run_type == "model-vs-model":
        ref_data = Dataset(parameter, data_type="ref")
        ref_use_cosp_variables = False

        try:
            ref_landfrac_ds = ref_data.get_time_series_dataset(
                "LANDFRAC", single_point=True
            )
            landfrac = ref_landfrac_ds["LANDFRAC"].sel(lat=slice(-70, -30))
        except Exception as e:
            logger.warning(f"Could not load reference LANDFRAC: {e}")
            landfrac = None

        try:
            # Try CLD_CAL_TMPICE first, fall back to CLDICE
            try:
                ref_cice_ds = ref_data.get_time_series_dataset(
                    "CLD_CAL_TMPICE", single_point=True
                )
                cice = ref_cice_ds["CLD_CAL_TMPICE"].sel(lat=slice(-70, -30))
                ref_use_cosp_variables = True
                logger.info("Using CLD_CAL_TMPICE for reference cloud ice")
            except Exception:
                ref_cice_ds = ref_data.get_time_series_dataset(
                    "CLDICE", single_point=True
                )
                cice = ref_cice_ds["CLDICE"].sel(lat=slice(-70, -30))
                logger.info("Using CLDICE for reference cloud ice")

            # Try CLD_CAL_TMPLIQ first, fall back to CLDLIQ
            try:
                ref_cliq_ds = ref_data.get_time_series_dataset(
                    "CLD_CAL_TMPLIQ", single_point=True
                )
                cliq = ref_cliq_ds["CLD_CAL_TMPLIQ"].sel(lat=slice(-70, -30))
                ref_use_cosp_variables = True
                logger.info("Using CLD_CAL_TMPLIQ for reference cloud liquid")
            except Exception:
                ref_cliq_ds = ref_data.get_time_series_dataset(
                    "CLDLIQ", single_point=True
                )
                cliq = ref_cliq_ds["CLDLIQ"].sel(lat=slice(-70, -30))
                logger.info("Using CLDLIQ for reference cloud liquid")

            # Load temperature: cosp_temp for COSP variables, T for original variables
            if ref_use_cosp_variables:
                try:
                    ref_temp_ds = ref_data.get_time_series_dataset(
                        "cosp_temp", single_point=True
                    )
                    temp = ref_temp_ds["cosp_temp"].sel(lat=slice(-70, -30))
                    logger.info("Using cosp_temp for reference temperature")
                except Exception:
                    logger.warning(
                        "Could not load reference cosp_temp, falling back to T"
                    )
                    ref_temp_ds = ref_data.get_time_series_dataset(
                        "T", single_point=True
                    )
                    temp = ref_temp_ds["T"].sel(lat=slice(-70, -30))
            else:
                ref_temp_ds = ref_data.get_time_series_dataset("T", single_point=True)
                temp = ref_temp_ds["T"].sel(lat=slice(-70, -30))
                logger.info("Using T for reference temperature")

        except Exception as e:
            logger.info(f"Error loading variables for reference data: {e}")
            raise

        parameter.ref_name_yrs = ref_data.get_name_yrs_attr(season)
        metrics_dict["ref"] = {}
        metrics_dict["ref"]["T"], metrics_dict["ref"]["LCF"] = compute_lcf(
            cice, cliq, temp, landfrac
        )
    parameter.output_file = "mixed-phase_partition"

    # TODO: save metrics
    plot(metrics_dict, parameter)

    return parameter
