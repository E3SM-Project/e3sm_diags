import os
import subprocess
import sys
import traceback
from datetime import datetime
from typing import TypedDict

import dask
import dask.bag as db

import e3sm_diags
from e3sm_diags.logger import LOG_FILENAME, _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter
from e3sm_diags.parser import SET_TO_PARSER
from e3sm_diags.parser.core_parser import CoreParser
from e3sm_diags.viewer.main import create_viewer

logger = _setup_child_logger(__name__)


class ProvPaths(TypedDict):
    """
    ProvPaths is a TypedDict that defines the structure for provenance paths.

    Attributes
    ----------
    results_dir: str
        Path to the diagnostic results.
    log_path : str
        Path to the log directory.
    parameter_files_path : str
        Path to the parameter files.
    python_script_path : str
        Path to the Python script.
    env_yml_path : str
        Path to the environment YAML file.
    index_html_path : str
        Path to the provenance index HTML file.
    """

    results_dir: str
    log_path: str
    parameter_files_path: str | None
    python_script_path: str | None
    env_yml_path: str | None
    index_html_path: str | None


def get_default_diags_path(set_name, run_type, print_path=True):
    """
    Returns the path for the default diags for plotset set_name.
    These are different depending on the run_type.
    """
    folder = "{}".format(set_name)
    fnm = "{}_{}.cfg".format(set_name, run_type)
    pth = os.path.join(e3sm_diags.INSTALL_PATH, folder, fnm)

    if print_path:
        logger.info("Using {} for {}.".format(pth, set_name))
    if not os.path.exists(pth):
        raise RuntimeError(
            "Plotting via set '{}' not supported, file {} not installed".format(
                set_name, fnm
            )
        )
    return pth


def save_provenance(results_dir: str, parser: CoreParser, no_viewer: bool) -> ProvPaths:
    """
    Store the provenance in results_dir.
    """
    prov_dir = os.path.join(results_dir, "prov")

    paths: ProvPaths = {
        "results_dir": results_dir,
        "log_path": os.path.join(prov_dir, LOG_FILENAME),
        "parameter_files_path": None,
        "python_script_path": None,
        "env_yml_path": None,
        "index_html_path": None,
    }

    if no_viewer:
        paths.update(
            {
                "parameter_files_path": "N/A (No viewer generated)",
                "python_script_path": "N/A (No viewer generated)",
                "env_yml_path": "N/A (No viewer generated)",
                "index_html_path": "N/A (No viewer generated)",
            }
        )

        return paths

    paths["parameter_files_path"] = _save_parameter_files(prov_dir, parser)
    paths["python_script_path"] = _save_python_script(prov_dir, parser)

    # FIXME: Replace Exception with specific exception type.
    try:
        paths["env_yml_path"] = _save_env_yml(prov_dir)
    except Exception:
        paths["env_yml_path"] = None
        traceback.print_exc()

    if not os.path.exists(prov_dir):
        os.makedirs(prov_dir, 0o755)

    # Create an HTML file to list the contents of the prov dir.
    index_html_path = os.path.join(prov_dir, "index.html")
    paths["index_html_path"] = index_html_path

    with open(index_html_path, "w") as f:
        f.write("<html><body><h1>Provenance Files</h1><ul>")

        for root, _, files in os.walk(prov_dir):
            for file_name in files:
                file_path = os.path.relpath(os.path.join(root, file_name), prov_dir)
                f.write(
                    f'<li><a href="{file_path}" target="_blank">{file_name}</a></li>'
                )

        f.write("</ul></body></html>")

    return paths


def _save_env_yml(results_dir: str) -> str | None:
    """
    Save the yml to recreate the environment in results_dir.
    """
    cmd = "conda env export"
    p = subprocess.Popen(cmd.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, err = p.communicate()

    filename = None

    if err:
        logger.exception("Error when creating env yml file: ")
        logger.exception(err)
    else:
        filename = os.path.join(results_dir, "environment.yml")

        with open(filename, "w") as f:
            f.write(output.decode("utf-8"))

    return filename


def _save_parameter_files(results_dir: str, parser: CoreParser) -> str | None:
    """
    Save the command line arguments used, and any py or cfg files.
    """
    filepath = os.path.join(results_dir, "cmd_used.txt")
    new_filepath = None

    cmd_used = " ".join(sys.argv)
    with open(filepath, "w") as f:
        f.write(cmd_used)

    args = parser.view_args()

    if hasattr(args, "parameters") and args.parameters:
        filepath = args.parameters
    elif hasattr(args, "other_parameters") and args.other_parameters:
        filepath = args.other_parameters[0]

    if not os.path.isfile(filepath):
        logger.warning("File does not exist: {}".format(filepath))
    else:
        with open(filepath, "r") as f:
            contents = "".join(f.readlines())

        # Remove any path, just keep the filename.
        new_filepath = filepath.split("/")[-1]
        new_filepath = os.path.join(results_dir, new_filepath)

        with open(new_filepath, "w") as f:
            f.write(contents)

    return new_filepath


def _save_python_script(results_dir: str, parser: CoreParser) -> str | None:
    """
    When using a Python script to run the
    diags via the API, dump a copy of the script.
    """
    args = parser.view_args()

    # FIXME: Is this code still needed?
    # If running the legacy way, there's nothing to be saved.
    if args.parameters:
        return None

    # Get the last argument that has .py in it.
    py_files = [f for f in sys.argv if f.endswith(".py")]

    # User didn't pass in a Python file, so they maybe ran:
    #    e3sm_diags -d diags.cfg
    if not py_files:
        return None

    fnm = py_files[-1]

    if not os.path.isfile(fnm):
        logger.warning("File does not exist: {}".format(fnm))
        return None

    with open(fnm, "r") as f:
        contents = "".join(f.readlines())

    # Remove any path, just keep the filename.
    new_filepath = fnm.split("/")[-1]
    new_filepath = os.path.join(results_dir, new_filepath)

    with open(new_filepath, "w") as f:
        f.write(contents)

    return new_filepath


# FIXME: B008 Do not perform function call `CoreParser` in argument defaults;
# instead, perform the call within the function, or read the default from a
# module-level singleton variable
def get_parameters(parser=CoreParser()):  # noqa B008
    """
    Get the parameters from the parser.
    """
    # A separate parser to just get the args used.
    # The reason it's a separate object than `parser`
    # is so we can parse the known args.
    parser_for_args = CoreParser()
    # The unknown args are _.
    # These are any set-specific args that aren't needed
    # for now, we just want to know what args are used.
    args, _ = parser_for_args.parser.parse_known_args()

    # Below is the legacy way to run this software, pre v2.0.0.
    # There weren't any arguments defined.
    if not any(getattr(args, arg) for arg in vars(args)):
        parser.print_help()
        sys.exit()

    # For when a user runs the software with commands like:
    #    e3sm_diags lat_lon [the other parameters]
    # This use-case is usually ran when the provenance
    # command is copied and pasted from the viewers.
    if args.set_name in SET_TO_PARSER:
        parser = SET_TO_PARSER[args.set_name]()
        parameters = parser.get_parameters(
            cmd_default_vars=False, argparse_vals_only=False
        )

    # The below two clauses are for the legacy way to
    # run this software, pre v2.0.0.
    # Ex: e3sm_diags -p params.py -d diags.cfg
    elif args.parameters and not args.other_parameters:  # -p only
        original_parameter = parser.get_orig_parameters(argparse_vals_only=False)

        # Load the default cfg files.
        run_type = getattr(original_parameter, "run_type", "model_vs_obs")
        default_diags_paths = [
            get_default_diags_path(set_name, run_type)
            for set_name in CoreParameter().sets
        ]

        other_parameters = parser.get_cfg_parameters(
            files_to_open=default_diags_paths, argparse_vals_only=False
        )

        parameters = parser.get_parameters(
            orig_parameters=original_parameter,
            other_parameters=other_parameters,
            cmd_default_vars=False,
            argparse_vals_only=False,
        )

    else:
        parameters = parser.get_parameters(
            cmd_default_vars=False, argparse_vals_only=False
        )

    parser.check_values_of_params(parameters)

    if not parameters:
        msg = "No parameters were able to be created. Please check your .py "
        msg += "file, and any .cfg files or command line args you're using."
        raise RuntimeError(msg)

    return parameters


def create_parameter_dict(parameters):
    d: dict[type, int] = dict()
    for parameter in parameters:
        t = type(parameter)
        if t in d.keys():
            d[t] += 1
        else:
            d[t] = 1
    return d


def _run_serially(parameters: list[CoreParameter]) -> list[CoreParameter]:
    """Run diagnostics with the parameters serially.

    Parameters
    ----------
    parameters : list[CoreParameter]
        The list of CoreParameter objects to run diagnostics on.

    Returns
    -------
    list[CoreParameter]
        The list of CoreParameter objects with results from the diagnostic run.
    """
    # A nested list of lists, where a sub-list represents the results of
    # the sets related to the CoreParameter object.
    nested_results: list[list[CoreParameter]] = []

    for parameter in parameters:
        nested_results.append(parameter._run_diag())

    # `results` becomes a list of lists of parameters so it needs to be
    # collapsed a level.
    collapsed_results = _collapse_results(nested_results)

    return collapsed_results


def _run_with_dask(parameters: list[CoreParameter]) -> list[CoreParameter]:
    """Run diagnostics with the parameters in parallel using Dask.

    This function passes ``run_diag`` to ``dask.bag.map``, which gets executed
    in parallel with ``.compute``.

    The first CoreParameter object's `num_workers` attribute is used to set
    the number of workers for ``.compute``.

    Parameters
    ----------
    parameters : list[CoreParameter]
        The list of CoreParameter objects to run diagnostics on.

    Returns
    -------
    list[CoreParameter]
        The list of CoreParameter objects with results from the diagnostic run.

    Notes
    -----
    https://docs.dask.org/en/stable/generated/dask.bag.map.html
    https://docs.dask.org/en/stable/generated/dask.dataframe.DataFrame.compute.html
    """
    bag = db.from_sequence(parameters)
    config = {"scheduler": "processes", "multiprocessing.context": "fork"}

    num_workers = getattr(parameters[0], "num_workers", None)
    if num_workers is None:
        raise ValueError(
            "The `num_workers` attribute is required for multiprocessing but it is not "
            "defined on the CoreParameter object. Set this attribute and try running "
            "again."
        )

    with dask.config.set(config):
        results = bag.map(CoreParameter._run_diag).compute(num_workers=num_workers)

    # `results` becomes a list of lists of parameters so it needs to be
    # collapsed a level.
    collapsed_results = _collapse_results(results)

    return collapsed_results


def _collapse_results(parameters: list[list[CoreParameter]]) -> list[CoreParameter]:
    """Collapses the results of diagnostic runs by one list level.

    Parameters
    ----------
    parameters : list[list[CoreParameter]]
        A list of lists of CoreParameter objects with results from the
        diagnostic run.

    Returns
    -------
    list[CoreParameter]
        A list of CoreParameter objects with results from the diagnostic run.
    """
    output_parameters = []

    for p1 in parameters:
        if isinstance(p1, list):
            for p2 in p1:
                output_parameters.append(p2)
        else:
            output_parameters.append(p1)

    return output_parameters


# FIXME: B006 Do not use mutable data structures for argument defaults
def main(parameters=[]) -> list[CoreParameter]:  # noqa B006
    # Get the diagnostic run parameters
    # ---------------------------------
    parser = CoreParser()

    # If no parameters are passed, use the parser args as defaults. Otherwise,
    # create the dictionary of expected parameters.
    if len(parameters) == 0:
        parameters = get_parameters(parser)

    expected_parameters = create_parameter_dict(parameters)

    if not os.path.exists(parameters[0].results_dir):
        os.makedirs(parameters[0].results_dir, 0o755)

    prov_paths = save_provenance(
        parameters[0].results_dir, parser, parameters[0].no_viewer
    )
    _log_diagnostic_run_info(prov_paths)

    # Perform the diagnostic run
    # --------------------------
    if parameters[0].multiprocessing:
        parameters_results = _run_with_dask(parameters)
    else:
        parameters_results = _run_serially(parameters)

    # Generate the viewer outputs using results
    # -----------------------------------------
    if not parameters_results:
        logger.warning(
            "There was not a single valid diagnostics run, no viewer created."
        )
    else:
        # If you get `AttributeError: 'NoneType' object has no attribute
        # `'no_viewer'` on this line then `run_diag` likely returns `None`.
        if parameters_results[0].no_viewer:
            logger.info("Viewer not created because the no_viewer parameter is True.")
        else:
            index_path = create_viewer(parameters_results)
            logger.info("Viewer HTML generated at {}".format(index_path))

    # Validate actual and expected parameters are aligned
    # ---------------------------------------------------
    actual_parameters = create_parameter_dict(parameters_results)
    if parameters_results[0].fail_on_incomplete and (
        actual_parameters != expected_parameters
    ):
        d: dict[type, tuple[int, int]] = dict()

        # Loop through all expected parameter types.
        for t in expected_parameters.keys():
            d[t] = (actual_parameters[t], expected_parameters[t])

        message = (
            "Not all parameters completed successfully. Check output above for "
            "errors/exceptions. The following dictionary maps parameter types to their "
            f"actual and expected numbers: {d}"
        )
        raise Exception(message)

    return parameters_results


def _log_diagnostic_run_info(prov_paths: ProvPaths):
    """Logs information about the diagnostic run.

    This method is useful for tracking the provenance of the diagnostic run
    and understanding the context of the diagnostic results.

    It logs the following information:
        - Timestamp of the run
        - Version information (Git branch and commit hash or module version)
        - Paths to the provenance files (log, parameter files, Python script,
          env yml, index HTML)

    Parameters
    ----------
    prov_paths : ProvPaths
        The paths to the provenance files.

    Notes
    -----
    The version information is retrieved from the current Git branch and
    commit hash. If the Git information is not available, it falls back
    to the version defined in the `e3sm_diags` module.
    """
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    try:
        branch_name = (
            subprocess.check_output(
                ["git", "rev-parse", "--abbrev-ref", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=subprocess.DEVNULL,
            )
            .strip()
            .decode("utf-8")
        )
        commit_hash = (
            subprocess.check_output(
                ["git", "rev-parse", "HEAD"],
                cwd=os.path.dirname(__file__),
                stderr=subprocess.DEVNULL,
            )
            .strip()
            .decode("utf-8")
        )
        version_info = f"branch {branch_name} with commit {commit_hash}"
    except subprocess.CalledProcessError:
        version_info = f"version {e3sm_diags.__version__}"

    (
        results_dir,
        log_path,
        parameter_files_path,
        python_script_path,
        env_yml_path,
        index_html_path,
    ) = prov_paths.values()

    logger.info(
        f"\n{'=' * 80}\n"
        f"E3SM Diagnostics Run\n"
        f"{'-' * 20}\n"
        f"Timestamp: {timestamp}\n"
        f"Version Info: {version_info}\n"
        f"Results Path: {results_dir}\n"
        f"Log Path: {log_path}\n"
        f"Parameter Files Path: {parameter_files_path}\n"
        f"Python Script Path: {python_script_path}\n"
        f"Environment YML Path: {env_yml_path}\n"
        f"Provenance Index HTML Path: {index_html_path}\n"
        f"{'=' * 80}\n"
    )


if __name__ == "__main__":
    main()
