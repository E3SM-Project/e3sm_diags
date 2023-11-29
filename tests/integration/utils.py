import ast
import configparser
import os
import subprocess
from typing import List

from e3sm_diags.parameter import SET_TO_PARAMETERS
from e3sm_diags.parameter.area_mean_time_series_parameter import (
    AreaMeanTimeSeriesParameter,
)
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parameter.enso_diags_parameter import EnsoDiagsParameter
from e3sm_diags.parameter.meridional_mean_2d_parameter import MeridionalMean2dParameter
from e3sm_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter


def run_cmd_and_pipe_stderr(command: str) -> List[str]:
    """Runs the test command and pipes the stderr for further processing.

    E3SM diags uses the Python logging module for logging runs. The Python
    logger uses stderr for streaming, rather than stdout (e.g., print
    statements). To capture information such as file paths (e.g.,
    ``reference_dir``) from the logger, stderr must be piped using
    ``capture_output=True``. Be aware, piping stderr results in logger messages
    not outputting to the console when running tests. The workaround is to
    perform a normal print of the entire piped stderr outputs once testing
    complete.

    Parameters
    ----------
    command : str
        The test command.

    Returns
    -------
    List[str]
        List of strings from stderr, decoded with "utf-8".

    Notes
    -----
    If capture_output is true, stdout and stderr will be captured. When used,
    the internal Popen object is automatically created with stdout=PIPE and
    stderr=PIPE. The stdout and stderr arguments may not be supplied at the
    same time as capture_output. If you wish to capture and combine both streams
    into one, use stdout=PIPE and stderr=STDOUT instead of capture_output.

    References
    ----------
    https://docs.python.org/3/library/subprocess.html
    """
    print("\nRunning tests, please wait for log output.")
    proc: subprocess.CompletedProcess = subprocess.run(
        command.split(), capture_output=True
    )
    stderr = proc.stderr.decode("utf-8").splitlines()

    print(*stderr, sep="\n")
    return stderr


def _get_test_params() -> List[CoreParameter]:
    param = CoreParameter()
    ts_param = AreaMeanTimeSeriesParameter()

    m2d_param = MeridionalMean2dParameter()
    m2d_param.plevs = [
        200.0,
        500.0,
    ]
    z2d_param = ZonalMean2dParameter()
    z2d_param.plevs = [
        200.0,
        300.0,
    ]

    enso_param = EnsoDiagsParameter()
    enso_param.test_name = "e3sm_v1"

    params = [param, ts_param, m2d_param, z2d_param, enso_param]

    return params


def _convert_cfg_to_param_objs(cfg_path: str) -> List[CoreParameter]:
    """Convert diagnostic cfg entries to parameter objects.

    NOTE: ast.literal_eval is not considered "safe" on untrusted data.
    The reason why it is used is because `configparser.ConfigParser`
    doesn't work well with parsing Python types from strings in
    `.cfg` files, resulting in things such as nested strings or string
    representation of lists. Since we are only calling literal_eval on
    `.cfg` files hosted in this repo, there is minimal risk here.

    Returns
    -------
    List[CoreParameter]
        A list of CoreParameter objects, one for each diagnotic set.
    """
    config = configparser.ConfigParser()
    config.read(cfg_path)
    params = []

    for set_name in config.sections():
        param = SET_TO_PARAMETERS[set_name]()

        for option in config.options(set_name):
            val = config.get(set_name, option)
            val = ast.literal_eval(val)

            setattr(param, option, val)

        params.append(param)

    return params


def _count_images(directory: str):
    """Count the number of images of type file_type in directory"""
    count = 0

    for _, __, files in os.walk(directory):
        for f in files:
            if f.endswith("png"):
                count += 1

    return count
