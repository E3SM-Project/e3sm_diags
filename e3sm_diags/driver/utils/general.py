import errno
import os

from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def strictly_increasing(L):
    return all(x < y for x, y in zip(L, L[1:]))


def strictly_decreasing(L):
    return all(x > y for x, y in zip(L, L[1:]))


def monotonically_decreasing(L):
    return all(x >= y for x, y in zip(L, L[1:]))


def monotonically_increasing(L):
    return all(x <= y for x, y in zip(L, L[1:]))


def monotonic(L):
    return monotonically_increasing(L) or monotonically_decreasing(L)


def get_output_dir(set_num, parameter):
    """
    Get the directory of where to save the outputs for a run.
    """
    results_dir = parameter.results_dir
    pth = os.path.join(results_dir, "{}".format(set_num), parameter.case_id)

    if not os.path.exists(pth):
        # When running diags in parallel, sometimes another process will create the dir.
        try:
            os.makedirs(pth, 0o755)
        except OSError as e:
            if e.errno != errno.EEXIST:
                raise

    return pth
