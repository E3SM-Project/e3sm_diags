from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def monotonic(L):
    return _monotonically_increasing(L) or _monotonically_decreasing(L)


def _monotonically_decreasing(L):
    return all(x >= y for x, y in zip(L, L[1:], strict=False))


def _monotonically_increasing(L):
    return all(x <= y for x, y in zip(L, L[1:], strict=False))
