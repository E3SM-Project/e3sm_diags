import pytest

from e3sm_diags.parser.area_mean_time_series_parser import AreaMeanTimeSeriesParser
from e3sm_diags.parser.arm_diags_parser import ARMDiagsParser
from e3sm_diags.parser.core_parser import CoreParser
from e3sm_diags.parser.diurnal_cycle_parser import DiurnalCycleParser
from e3sm_diags.parser.enso_diags_parser import EnsoDiagsParser
from e3sm_diags.parser.meridional_mean_2d_parser import MeridionalMean2dParser
from e3sm_diags.parser.qbo_parser import QboParser
from e3sm_diags.parser.streamflow_parser import StreamflowParser
from e3sm_diags.parser.tc_analysis_parser import TCAnalysisParser
from e3sm_diags.parser.zonal_mean_2d_parser import ZonalMean2dParser
from e3sm_diags.parser.zonal_mean_2d_stratosphere_parser import (
    ZonalMean2dStratosphereParser,
)


class TestCoreParser:
    def test__init__(self):
        CoreParser()

    @pytest.mark.xfail
    def test_check_values_of_params(self):
        assert 0

    @pytest.mark.xfail
    def test_add_arguments(self):
        assert 0

    @pytest.mark.xfail
    def test_parse_args(self):
        assert 0

    @pytest.mark.xfail
    def test_parse_args_removes_ipykernel_args(self):
        assert 0

    @pytest.mark.xfail
    def test_view_args_returns_parser_namespace(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_cmdline_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_orig_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_other_parameters(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_default_vars(self):
        assert 0

    @pytest.mark.xfail
    def test_get_parameters_returns_cmd_default_vars(self):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_returns_parameters_created_by_running_from_CLI(
        self,
    ):
        # FIXME: Should we deprecate this method for running `e3sm_diags`?
        # https://e3sm-project.github.io/e3sm_diags/_build/html/main/config-run.html#e3sm-diags-p-older-method
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_returns_parameters_created_by_cfg_file(self):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_checks_values_in_cfg_file_and_returns_parameter(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_get_cfg_parameters_only_uses_argparse_values_and_returns_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_cmdline_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_orig_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0

    @pytest.mark.xfail
    def test_select_returns_other_parameters_that_are_subset_of_the_main_parameters(
        self,
    ):
        assert 0


def test_area_mean_time_series_parser_initializes():
    AreaMeanTimeSeriesParser()


def test_arms_diags_parser_initializes():
    ARMDiagsParser()


def test_diurnal_cycle_parser_initializes():
    DiurnalCycleParser()


def test_enso_diags_parser_initializes():
    EnsoDiagsParser()


def test_meridional_mean_2d_parser_initializes():
    MeridionalMean2dParser()


def test_qbo_parser_initializes():
    QboParser()


def test_streamflow_parser_initializes():
    StreamflowParser()


def test_tc_analysis_parser_initializes():
    TCAnalysisParser()


def test_zonal_mean_2d_parser_initializes():
    ZonalMean2dParser()


@pytest.mark.xfail
def test_zonal_mean_2d_stratosphere_parser_initializes():
    # FAILED tests/e3sm_diags/test_parsers.py::test_zonal_mean_2d_stratosphere_parser_initializes - TypeError: e3sm_diags.parser.core_parser.CoreParser.__init__() got multiple values for keyword argument 'parameter_cls'
    # FIXME: This class has multiple inheritance (ZonalMean2dParser -> CoreParser)
    ZonalMean2dStratosphereParser()
