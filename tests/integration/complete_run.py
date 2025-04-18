"""The complete E3SM Diagnostic run.

Due to the large amount of data required to run, this test will be run manually
on Anvil (rather than as part of the CI tests).

This test should be run with the latest E3SM Diags tutorial code.

Run the following first:
  - srun --pty --nodes=1 --time=01:00:00 /bin/bash
  - source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
  - Or: source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh
"""

import os

from e3sm_diags.parameter.annual_cycle_zonal_mean_parameter import ACzonalmeanParameter
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.arm_diags_parameter import ARMDiagsParameter
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.diurnal_cycle_parameter import DiurnalCycleParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.qbo_parameter import QboParameter
from e3sm_diags.parameter.streamflow_parameter import StreamflowParameter
from e3sm_diags.parameter.tc_analysis_parameter import TCAnalysisParameter
from e3sm_diags.parameter.zonal_mean_2d_stratosphere_parameter import (
    ZonalMean2dStratosphereParameter,
)
from e3sm_diags.run import runner
from tests.integration.utils import _compare_images


class TestCompleteRun:
    def test_complete_run(self):
        actual_images_dir = run_lcrc(".")

        # The expected_images_file lists all images we expect to compare.
        expected_images_file = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/image_list_all_sets.txt"
        expected_images_dir = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/unit_test_complete_run/expected/all_sets"

        mismatched_images = []  # type:ignore

        with open(expected_images_file) as f:
            for line in f:
                image_name = line.strip("./").strip("\n")
                path_to_actual_png = os.path.join(actual_images_dir, image_name)
                path_to_expected_png = os.path.join(expected_images_dir, image_name)

                mismatched_images = _compare_images(
                    mismatched_images,
                    image_name,
                    path_to_actual_png,
                    path_to_expected_png,
                )

        assert len(mismatched_images) == 0


def run_lcrc(html_prefix):
    # Run the following first:
    # srun --pty --nodes=1 --time=01:00:00 /bin/bash
    # source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
    # Or: source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh

    ref_data_prefix = "/lcrc/group/e3sm/public_html/diagnostics/observations/Atm"
    test_data_prefix = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
    # use highrequency grid box output at ARM sites from another simulation when the output is available
    test_data_prefix2 = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags/20210719.PhaseII.F20TR-P3.NGD.ne30pg2.compy"

    d = dict()

    d["obs_climo"] = os.path.join(ref_data_prefix, "climatology/")
    d["test_climo"] = os.path.join(test_data_prefix, "climatology/rgr/")

    d["obs_ts"] = os.path.join(ref_data_prefix, "time-series/")
    d["test_ts"] = os.path.join(test_data_prefix, "time-series/rgr/")

    d["dc_obs_climo"] = os.path.join(ref_data_prefix, "climatology/")
    d["dc_test_climo"] = os.path.join(test_data_prefix, "diurnal_climatology/rgr")

    d["arm_obs"] = os.path.join(ref_data_prefix, "arm-diags-data/")
    d["arm_test"] = os.path.join(test_data_prefix2, "arm-diags-data/")

    d["tc_obs"] = os.path.join(ref_data_prefix, "tc-analysis/")
    d["tc_test"] = os.path.join(test_data_prefix, "tc-analysis/")

    return run_all_sets(html_prefix, d)


def run_all_sets(html_prefix, d):
    param = CoreParameter()

    param.reference_data_path = d["obs_climo"]
    param.test_data_path = d["test_climo"]
    param.test_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
    param.seasons = [
        "ANN",
        "JJA",
    ]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

    param.results_dir = os.path.join(html_prefix, "v2_6_0_all_sets")
    param.multiprocessing = True
    param.num_workers = 5

    # Set specific parameters for new sets
    enso_param = EnsoDiagsParameter()
    enso_param.reference_data_path = d["obs_ts"]
    enso_param.test_data_path = d["test_ts"]
    enso_param.test_name = "e3sm_v2"
    enso_param.test_start_yr = "0051"
    enso_param.test_end_yr = "0060"
    # Enso obs data range from year 1979 to 2016
    enso_param.ref_start_yr = "2001"
    enso_param.ref_end_yr = "2010"

    qbo_param = QboParameter()
    qbo_param.reference_data_path = d["obs_ts"]
    qbo_param.test_data_path = d["test_ts"]
    qbo_param.test_name = "e3sm_v2"
    qbo_param.start_yr = "0051"
    qbo_param.end_yr = "0060"
    # Qbo obs data range from year 1979 to 2019
    # Number of years of test and ref should match
    qbo_param.ref_start_yr = "2001"
    qbo_param.ref_end_yr = "2010"

    ts_param = AreaMeanTimeSeriesParameter()
    ts_param.reference_data_path = d["obs_ts"]
    ts_param.test_data_path = d["test_ts"]
    ts_param.test_name = "e3sm_v2"
    ts_param.start_yr = "0051"
    ts_param.end_yr = "0060"

    dc_param = DiurnalCycleParameter()
    dc_param.reference_data_path = d["dc_obs_climo"]
    dc_param.test_data_path = d["dc_test_climo"]
    dc_param.short_test_name = "e3sm_v2"
    # Plotting diurnal cycle amplitude on different scales. Default is True
    dc_param.normalize_test_amp = False

    streamflow_param = StreamflowParameter()
    streamflow_param.reference_data_path = d["obs_ts"]
    streamflow_param.test_data_path = d["test_ts"]
    streamflow_param.short_test_name = "e3sm_v2"
    streamflow_param.test_start_yr = "0051"
    streamflow_param.test_end_yr = "0060"
    # Streamflow gauge station data range from year 1986 to 1995
    streamflow_param.ref_start_yr = "1986"
    streamflow_param.ref_end_yr = "1995"

    arm_param = ARMDiagsParameter()
    arm_param.reference_data_path = d["arm_obs"]
    arm_param.ref_name = "armdiags"
    arm_param.test_data_path = d["arm_test"]
    arm_param.test_name = "e3sm_v2"
    arm_param.test_start_yr = "1996"
    arm_param.test_end_yr = "2010"
    # For model vs obs, the ref start and end year can be any four digit strings for now, will use all available years form obs
    arm_param.ref_start_yr = "0001"
    arm_param.ref_end_yr = "0001"

    tc_param = TCAnalysisParameter()
    tc_param.reference_data_path = d["tc_obs"]
    tc_param.test_data_path = d["tc_test"]
    tc_param.short_test_name = "e3sm_v2"
    tc_param.test_start_yr = "0051"
    tc_param.test_end_yr = "0060"
    # For model vs obs, the ref start and end year can be any four digit strings for now, use all available years form obs by default
    tc_param.ref_start_yr = "1979"
    tc_param.ref_end_yr = "2018"

    ac_param = ACzonalmeanParameter()
    zm_param = ZonalMean2dStratosphereParameter()

    runner.sets_to_run = [
        "lat_lon",
        "zonal_mean_xy",
        "zonal_mean_2d",
        "zonal_mean_2d_stratosphere",
        "polar",
        "cosp_histogram",
        "meridional_mean_2d",
        "annual_cycle_zonal_mean",
        "enso_diags",
        "qbo",
        "area_mean_time_series",
        "diurnal_cycle",
        "streamflow",
        "arm_diags",
        "tc_analysis",
    ]
    runner.run_diags(
        [
            param,
            zm_param,
            ac_param,
            enso_param,
            qbo_param,
            ts_param,
            dc_param,
            streamflow_param,
            arm_param,
            tc_param,
        ]
    )

    return param.results_dir
