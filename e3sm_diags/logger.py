"""Logger module for setting up a custom logger."""

import logging
import logging.handlers

LOG_FILENAME = "e3sm_diags_run.log"
LOG_FORMAT = (
    "%(asctime)s [%(levelname)s]: %(filename)s(%(funcName)s:%(lineno)s) >> %(message)s"
)
LOG_FILEMODE = "w"
LOG_LEVEL = logging.INFO

# Add a console handler to display warnings in the console. This is useful
# for when other package loggers raise warnings (e.g, NumPy, Xarray).
console_handler = logging.StreamHandler()
console_handler.setLevel(logging.INFO)
console_handler.setFormatter(logging.Formatter(LOG_FORMAT))
logging.getLogger().addHandler(console_handler)


def _setup_root_logger():
    """Configures the root logger.

    This function sets up the root logger with a predefined format and log level.
    It also enables capturing of warnings issued by the `warnings` module and
    redirects them to the logging system.

    Notes
    -----
    - The `force=True` parameter ensures that any existing logging configuration
      is overridden.
    - The file handler is added dynamically to the root logger later in the
      ``Run`` class once the log file path is known.
    """
    logging.basicConfig(
        format=LOG_FORMAT,
        level=LOG_LEVEL,
        force=True,
    )

    logging.captureWarnings(True)


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


def _add_filehandler(log_path: str):
    """Adds a file handler to the root logger dynamically.

    Adding the file handler will also create the log file automatically.

    Parameters
    ----------
    log_path : str
        The path to the log file.

    Notes
    -----
    Any warnings that appear before the log filehandler is instantiated will not
    be captured (e.g,. esmpy VersionWarning). However, they will still be
    captured by the console via the default StreamHandler.
    """
    file_handler = logging.FileHandler(log_path, mode=LOG_FILEMODE)
    file_handler.setFormatter(logging.Formatter(LOG_FORMAT))
    logging.root.addHandler(file_handler)
