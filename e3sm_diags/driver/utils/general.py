from e3sm_diags.logger import _setup_child_logger

logger = _setup_child_logger(__name__)


def monotonic(L):
    return _monotonically_increasing(L) or _monotonically_decreasing(L)


def _monotonically_decreasing(L):
    # FIXME: B905: zip() without an explicit strict= parameter
    return all(x >= y for x, y in zip(L, L[1:]))


def _monotonically_increasing(L):
    # FIXME: B905: zip() without an explicit strict= parameter
    return all(x <= y for x, y in zip(L, L[1:]))
