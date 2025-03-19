"""Logger module for setting up a custom logger."""

import logging
import logging.handlers
import os
import shutil

LOG_FILENAME = "e3sm_diags_run.log"
LOG_FORMAT = (
    "%(asctime)s [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >> %(message)s"
)
LOG_FILEMODE = "w"
LOG_LEVEL = logging.INFO


# Setup the root logger with a default log file.
# `force` is set to `True` to automatically remove root handlers whenever
# `basicConfig` called. This is required for cases where multiple e3sm_diags
# runs are executed. Otherwise, the logger objects attempt to share the same
# root file reference (which gets deleted between runs), resulting in
# `FileNotFoundError: [Errno 2] No such file or directory: 'e3sm_diags_run.log'`.
# More info here: https://stackoverflow.com/a/49202811
logging.basicConfig(
    format=LOG_FORMAT,
    filename=LOG_FILENAME,
    filemode=LOG_FILEMODE,
    level=LOG_LEVEL,
    force=True,
)
logging.captureWarnings(True)

# Add a console handler to display warnings in the console. This is useful
# for when other package loggers raise warnings (e.g, NumPy, Xarray).
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter(LOG_FORMAT))
logging.getLogger().addHandler(console_handler)


def _setup_child_logger(name: str, propagate: bool = True) -> logging.Logger:
    """Sets up a logger that is a child of the root logger.

    This child logger inherits the root logger's handlers.

    Parameters
    ----------
    name : str
        Name of the file where this function is called.
    propagate : bool, optional
        Whether to propagate logger messages or not, by default True.

    Returns
    -------
    logging.Logger
        The logger.

    Examples
    ---------
    Detailed information, typically of interest only when diagnosing problems:

    >>> logger.debug("")

    Confirmation that things are working as expected:

    >>> logger.info("")

    An indication that something unexpected happened, or indicative of some
    problem in the near future:

    >>> logger.warning("")

    The software has not been able to perform some function due to a more
    serious problem:

    >>> logger.error("")

    Similar to ``logger.error()``, but also outputs stack trace:

    >>> logger.exception("", exc_info=True)

    A serious error, indicating that the program itself may be unable to
    continue running:

    >>> logger.critical("")
    """
    logger = logging.getLogger(name)
    logger.propagate = propagate

    return logger


def _update_root_logger_filepath_to_prov_dir(log_path: str):
    """Updates the log file path to the provenance directory.

    This method changes the log file path to a subdirectory named 'prov'
    within the given results directory. It updates the filename of the
    existing file handler to the new path.

    Parameters
    ----------
    log_path : str
        The path to the log file, which is stored in the `results_dir`
        sub-directory called "prov".

    Notes
    -----
    - The method assumes that a logging file handler is already configured.
    - The log file is closed and reopened at the new location.
    - The log file mode is determined by the constant `LOG_FILEMODE`.
    - The log file name is determined by the constant `LOG_FILENAME`.
    """
    for handler in logging.root.handlers:
        if isinstance(handler, logging.FileHandler):
            handler.baseFilename = log_path
            handler.stream.close()
            handler.stream = open(log_path, LOG_FILEMODE)  # type: ignore
            break


def move_log_to_prov_dir(results_dir: str, logger: logging.Logger):
    """Moves the e3sm diags log file to the provenance directory.

    This function should be called at the end of the diagnostic run to capture
    all console outputs.

    Parameters
    ----------
    results_dir : str
        The results directory for the run.
    """
    provenance_dir = f"{results_dir}/prov/{LOG_FILENAME}"

    # Must copy and then delete because shutil.move does not work if different
    # filesystems are used for the source and destination directories.
    shutil.copy(LOG_FILENAME, provenance_dir)
    os.remove(LOG_FILENAME)
    logger.info(f"Log file saved in {provenance_dir}")
