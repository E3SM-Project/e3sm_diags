"""
Make sure to run the machine-specific commands below before
running this script:

Compy:
    srun --pty --nodes=1 --time=01:00:00 /bin/bash
    source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh

LCRC:
    srun --pty --nodes=1 --time=01:00:00 /bin/bash
    source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh
    Or: source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh

NERSC perlmutter cpu:
    salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm
    source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
"""
# flake8: noqa E501

from typing import Tuple, TypedDict

from mache import MachineInfo

from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.run import runner


class MachinePaths(TypedDict):
    html_path: str
    obs_climo: str
    test_climo: str
    obs_ts: str
    test_ts: str
    dc_obs_climo: str
    dc_test_climo: str
    arm_obs: str
    arm_test: str
    tc_obs: str
    tc_test: str


def run_all_sets():
    machine_paths: MachinePaths = _get_machine_paths()

    param = CoreParameter()

    param.reference_data_path = machine_paths["obs_climo"]
    param.test_data_path = machine_paths["test_climo"]
    param.test_name = "20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis"
    param.seasons = [
        "ANN",
        "JJA",
        "DJF",
    ]  # Default setting: seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

    param.results_dir = f"{machine_paths['html_path']}/v3_2_0_cosp_v2_login"
    param.multiprocessing = True
    param.num_workers = 24



    runner.sets_to_run = [
        "lat_lon",
        "cosp_histogram",
    ]

    runner.run_diags(
        [
            param,
        ]
    )

    return param.results_dir


def _get_machine_paths() -> MachinePaths:
    """Returns the paths on the machine that are required to run e3sm_diags.

    Returns
    -------
    MachinePaths
        A dictionary of paths on the machine, with the key being the path type
        and the value being the absolute path string.
    """
    # Get the current machine's configuration info.
    machine_info = MachineInfo()
    machine = machine_info.machine

    if machine not in ["anvil", "chrysalis", "compy", "pm-cpu", "cori-haswell", "cori-knl"]:
        raise ValueError(f"e3sm_diags is not supported on this machine ({machine}).")

    # Path to the HTML outputs for the current user.
    web_portal_base_path = machine_info.config.get("web_portal", "base_path")
    html_path = f"{web_portal_base_path}/{machine_info.username}/"

    # Path to the reference data directory.
    diags_base_path = machine_info.diagnostics_base
    ref_data_dir = f"{diags_base_path}/observations/Atm"

    # Paths to the test data directories.
    test_data_dir, test_data_dir2 = _get_test_data_dirs(machine)

    # Construct the paths required by e3sm_diags using the base paths above.
    machine_paths: MachinePaths = {
        "html_path": html_path,
        #"obs_climo": f"{ref_data_dir}/climatology",
        "obs_climo": "/pscratch/sd/c/chengzhu/tests/test_yuying_cosp",
        "test_climo": f"{test_data_dir}/climatology/rgr/",
        "obs_ts": f"{ref_data_dir}/time-series/",
        "test_ts": f"{test_data_dir}/time-series/rgr/",
        "dc_obs_climo": f"{ref_data_dir}/climatology",
        "dc_test_climo": f"{test_data_dir}/diurnal_climatology/rgr",
        "arm_obs": f"{ref_data_dir}/arm-diags-data/",
        "arm_test": f"{test_data_dir2}/arm-diags-data/",
        "tc_obs": f"{ref_data_dir}/tc-analysis/",
        "tc_test": f"{test_data_dir}/tc-analysis/",
    }

    return machine_paths


def _get_test_data_dirs(machine: str) -> Tuple[str, str]:
    """Get the directories for test data based on the machine.

    The second path is for using the high frequency grid box output at ARM sites
    from another simulation when the output is available.

    Parameters
    ----------
    machine : str
        The name of the machine.

    Returns
    -------
    Tuple[str, str]
        A tuple of two strings, each representing a test data directory path.
    """
    test_data_dirs = None

    # TODO: Update this function to use `mache` after the directories are updated.
    if machine in ["chrysalis", "anvil"]:
        base = "/lcrc/group/e3sm/public_html/e3sm_diags_test_data/postprocessed_e3sm_v2_data_for_e3sm_diags"
    elif machine in ["compy"]:
        base = "/compyfs/e3sm_diags_data/postprocessed_e3sm_v2_data_for_e3sm_diags"
    elif machine in ["cori-haswell", "cori-knl", "pm-cpu"]:
        base = "/global/cfs/cdirs/e3sm/e3sm_diags/postprocessed_e3sm_v2_data_for_e3sm_diags"

    test_data_dirs = (
        f"{base}/20210528.v2rc3e.piControl.ne30pg2_EC30to60E2r2.chrysalis",
        #f"{base}/20210719.PhaseII.F20TR-P3.NGD.ne30pg2.compy",
        f"{base}/20221103.v2.LR.amip.NGD_v3atm.chrysalis",
    )

    return test_data_dirs  # type: ignore


if __name__ == "__main__":
    run_all_sets()
