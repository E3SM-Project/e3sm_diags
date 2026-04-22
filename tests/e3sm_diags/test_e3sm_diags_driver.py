import pytest

from e3sm_diags.e3sm_diags_driver import _run_serially, _run_with_dask
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.parameter.core_parameter import CoreParameter

logger = _setup_child_logger("e3sm_diags.e3sm_diags_driver", propagate=True)


class TestRunDiag:
    @pytest.mark.xfail
    def test_run_diag_serially_returns_parameters_with_results(self):
        # FIXME: This test will fail while we refactor sets and utilities. It
        # should be fixed after all sets are refactored.
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = _run_serially([parameter])
        expected = [parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results == expected

    @pytest.mark.xfail
    def test_run_diag_with_dask_returns_parameters_with_results(self):
        # FIXME: This test will fail while we refactor sets and utilities. It
        # should be fixed after all sets are refactored.
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]

        results = _run_with_dask([parameter])

        expected_parameter = CoreParameter()
        expected_parameter.sets = ["lat_lon"]
        expected_parameter.current_set = "lat_lon"
        expected_parameter.test_name_yrs = ""
        expected_parameter.ref_name_yrs = ""
        expected_parameter.model_only = False
        expected = [expected_parameter]

        # NOTE: We are only testing that the function returns a list of
        # parameter objects, not the results themselves. There are integration
        # tests validates the results.
        assert results[0].__dict__ == expected[0].__dict__

    @pytest.mark.xfail
    def test_run_diag_with_dask_raises_error_if_num_workers_attr_not_set(
        self,
    ):
        # FIXME: This test will while we refactor sets and utilities. It should
        # be fixed after all sets are refactored.
        parameter = CoreParameter()
        parameter.sets = ["lat_lon"]
        del parameter.num_workers

        with pytest.raises(ValueError):
            _run_with_dask([parameter])
