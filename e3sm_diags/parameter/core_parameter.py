from __future__ import annotations

import copy
import importlib
import os
import re
import sys
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

from e3sm_diags.derivations.derivations import DerivedVariablesMap
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.general import pad_year
from e3sm_diags.driver.utils.regrid import REGRID_TOOLS
from e3sm_diags.logger import _setup_child_logger

if TYPE_CHECKING:
    import xarray as xr

    from e3sm_diags.driver.utils.type_annotations import TimeSelection, TimeSlice

logger = _setup_child_logger(__name__)

# FIXME: There is probably a better way of defining default sets because most of
# this is repeated in SETS_TO_PARAMETERS and SETS_TO_PARSERS.
# Also integration tests will break if "mp_partition" is included because
# we did not take it into account yet.
DEFAULT_SETS = [
    "zonal_mean_xy",
    "zonal_mean_2d",
    "zonal_mean_2d_stratosphere",
    "meridional_mean_2d",
    "lat_lon",
    "lat_lon_native",
    "polar",
    "area_mean_time_series",
    "cosp_histogram",
    "enso_diags",
    "qbo",
    "streamflow",
    "diurnal_cycle",
    "arm_diags",
    "tc_analysis",
    "annual_cycle_zonal_mean",
    "lat_lon_land",
    "lat_lon_river",
    "aerosol_aeronet",
    "aerosol_budget",
    #    "mp_partition",
    #    "tropical_subseasonal",
]

YEAR_ATTRIBUTES = [
    "start_yr",
    "end_yr",
    "test_start_yr",
    "test_end_yr",
    "ref_start_yr",
    "ref_end_yr",
]

if TYPE_CHECKING:
    from e3sm_diags.driver.utils.dataset_xr import Dataset


class CoreParameter:
    def __init__(self):
        # File I/O
        # ------------------------
        # FIXME: (REQUIRED) attributes should be set through `__init__` since
        # they are required, but they are currently set after initializing the
        # CoreParameter object (e.g., param.results_dir = "dir"). The default
        # values are set to "", and `self.check_values()` will validate that
        # it is not "" since they must be set by the user.

        # (REQUIRED) Path to the reference (obs) data. If there are multiple
        # datasets in the path, use ref_name parameter to specify which dataset
        # should be used.
        self.reference_data_path: str = ""

        # (REQUIRED) Path to the test (model) data.
        self.test_data_path: str = ""

        # File paths for data files (set dynamically during processing using
        # the dataset's file_path attribute).
        self.test_data_file_path = ""  # Path to test data file
        self.ref_data_file_path = ""  # Path to reference data file

        # (REQUIRED) The name of the folder where all runs will be stored.
        self.results_dir: str = ""

        # The name of the folder where the results (plots and nc files) will be
        # stored for a single run
        self.case_id: str = ""

        # Set to True to not generate a Viewer for the result.
        self.no_viewer: bool = False

        # If True, stops running all of the diagnostics on the first failure.
        # If False (the default), all errors are caught and ignored. If there
        # was an error and a plot could not be created, there’s a ‘—’ for that
        # set of parameters in the viewer.
        self.debug: bool = False

        # Time series data
        # ------------------------
        # Set to True if the ref data is in timeseries format. Default False. If
        # True, both ref_start_yr and ref_end_yr must also be set.
        self.ref_timeseries_input: bool = False

        # Set to True if the test data is in timeseries format. Default False.
        # If True, both test_start_yr and # test_end_yr must also be set.
        self.test_timeseries_input: bool = False

        # Diagnostic run settings
        # ------------------------
        # The supported run type for the diagnostics. Possible options are:
        # 'model_vs_obs' (by default), 'model_vs_model', or 'obs_vs_obs'.
        self.run_type: str = "model_vs_obs"

        # A list of the sets to be run, by default all sets.
        self.sets: list[str] = DEFAULT_SETS

        # The current set that is being ran when looping over sets in
        # `e3sm_diags_driver.run_diag()`.
        self.current_set: str = ""

        self.variables: list[str] = []
        self.seasons: list[ClimoFreq] = ["ANN", "DJF", "MAM", "JJA", "SON"]

        # Time slice parameters (mutually exclusive with seasons)
        # Index-based time selection for snapshot analysis using individual time indices
        # Examples: ["0"], ["5"], ["0", "1", "2"]
        self.time_slices: list[TimeSlice] = []

        self.regions: list[str] = ["global"]

        self.regrid_tool: REGRID_TOOLS = "xesmf"
        self.regrid_method: str = "conservative_normed"

        self.plevs: list[float] = []
        self.plot_log_plevs: bool = False
        self.plot_plevs: bool = False

        self.multiprocessing: bool = False
        self.num_workers: int = 4

        # Diagnostic plot settings
        # ------------------------
        self.main_title: str = ""
        # TODO: Remove `backend` because it is always e3sm_diags/plot/cartopy.
        # This change cascades down to changes in `e3sm_diags.plot.plot`.
        self.backend: str = "cartopy"
        self.save_netcdf: bool = False

        # Plot format settings
        self.output_file: str = ""
        self.output_format: list[str] = ["png"]
        self.output_format_subplot: list[str] = []
        self.canvas_size_w: int = 1212
        self.canvas_size_h: int = 1628
        self.figsize: tuple[float, float] = (8.5, 11.0)
        self.dpi: int = 150
        self.arrows: bool = True
        self.logo: bool = False
        self.contour_levels: list[float] = []

        # Test plot settings
        self.test_name: str = ""
        self.test_name_yrs: str = ""
        self.short_test_name: str = ""
        self.test_title: str = ""
        self.test_colormap: str = "cet_rainbow.rgb"
        self.test_units: str = ""

        # Reference plot settings
        # `ref_name` is used to search though the reference data directories.
        self.ref_name: str = ""
        self.ref_name_yrs: str = ""
        # `reference_name` is printed above ref plots.
        self.reference_name: str = ""
        self.short_ref_name: str = ""
        self.reference_title: str = ""
        self.reference_colormap: str = "cet_rainbow.rgb"
        self.reference_units: str = ""

        # The reference file name based on the season and other parameters, for
        # climo files only.
        self.ref_file: str = ""

        # Variable ID is used as the reference file ID in `general.save_ncfiles()`
        self.var_id: str = ""

        # Variable region, used in polar plots.
        self.var_region: str = ""

        # Difference plot settings
        self.diff_name: str = ""
        self.diff_title: str = "Model - Observation"
        self.diff_colormap: str = "diverging_bwr.rgb"
        self.diff_levels: list[float] = []
        self.diff_units: str = ""
        self.diff_type: str = "absolute"

        # Other settings
        # ------------------------
        # TODO: Need documentation on these attributes here and
        # here: https://e3sm-project.github.io/e3sm_diags/_build/html/main/available-parameters.html
        self.dataset: str = ""
        self.granulate: list[str] = [
            "variables",
            "seasons",
            "plevs",
            "regions",
            "time_slices",
        ]
        self.selectors: list[str] = ["sets", "seasons"]
        self.viewer_descr: dict[str, str] = {}
        self.fail_on_incomplete: bool = False

        # List of user derived variables, set in `dataset.Dataset`.
        self.derived_variables: DerivedVariablesMap = {}

        # FIXME: This attribute is only used in `lat_lon_driver.py`
        self.model_only: bool = False

    def __add__(self, other):
        """Copies attributes from the other CoreParameter object to this one.

        This method overrides existing attributes if they are already set
        using the new values from the other object.

        It is based on cdp.cdp_parameter.CDPParameter.__add__().

        Parameters
        ----------
        other : CoreParameter
            The other CoreParameter (or sub-class) object.

        Returns
        -------
        CoreParameter
            A new instance of CoreParameter with new attributes.
        """
        new_obj = copy.deepcopy(self)
        new_attrs = other.__dict__

        for attr, value in new_attrs.items():
            setattr(new_obj, attr, value)

        return new_obj

    def check_values(self):
        # The default values for these attributes are set to "" in `__init__`.
        # FIXME: The user should pass these values into the object
        # initialization rather than having to check they are set here.
        must_have_params = ["reference_data_path", "test_data_path", "results_dir"]
        for param in must_have_params:
            if getattr(self, param) == "":
                msg = "You need to specify {p} in the parameters file or via the command line using --{p}".format(
                    p=param
                )
                raise RuntimeError(msg)

        if self.ref_timeseries_input and not (
            hasattr(self, "ref_start_yr") and hasattr(self, "ref_end_yr")
        ):
            msg = (
                "You need to define both the 'ref_start_yr' and 'ref_end_yr' parameter."
            )
            raise RuntimeError(msg)

        if self.test_timeseries_input and not (
            hasattr(self, "test_start_yr") and hasattr(self, "test_end_yr")
        ):
            msg = "You need to define both the 'test_start_yr' and 'test_end_yr' parameter."
            raise RuntimeError(msg)

        if self.time_slices:
            self._validate_time_slice_format()

    def _validate_time_slice_format(self) -> None:
        """Validate that time_slice follows the expected format.

        Time slices must be non-negative integer indices representing
        individual time steps in the dataset.

        Parameters
        ----------
        time_slice : str
            The time slice string to validate. Must be a non-negative integer.

        Raises
        ------
        ValueError
            If the time slice format is invalid (not a non-negative integer).
        """
        # Define the regex pattern for a non-negative integer, including no
        # leading zeros except for zero itself.
        pattern = r"^(0|[1-9]\d*)$"

        for time_slice in self.time_slices:
            if not re.match(pattern, time_slice.strip()):
                raise ValueError(
                    f"Invalid time_slice format: '{time_slice}'. "
                    f"Expected a non-negative integer index. Examples: '0', '5', '42'"
                )

    def _set_param_output_attrs(
        self,
        var_key: str,
        time_selection: TimeSelection,
        region: str,
        ref_name: str,
        ilev: float | None,
    ):
        """Set the parameter output attributes based on argument values.

        Parameters
        ----------
        var_key : str
            The variable key.
        time_selection : TimeSelection
            The time slice or season.
        region : str
            The region.
        ref_name : str
            The reference name.
        ilev : float | None
            The pressure level, by default None. This option is only set if the
            variable is 3D.
        """
        if ilev is None:
            output_file = f"{ref_name}-{var_key}-{time_selection}-{region}"
            main_title = f"{var_key} {time_selection} {region}"
        else:
            ilev_str = str(int(ilev))
            output_file = f"{ref_name}-{var_key}-{ilev_str}-{time_selection}-{region}"
            main_title = f"{var_key} {ilev_str} mb {time_selection} {region}"

        self.output_file = output_file
        self.main_title = main_title

    def _set_name_yrs_attrs(
        self, ds_test: Dataset, ds_ref: Dataset, time_selection: TimeSelection | None
    ):
        """Set the test_name_yrs and ref_name_yrs attributes.

        Parameters
        ----------
        ds_test : Dataset
            The test dataset object used for setting ``self.test_name_yrs``.
        ds_ref : Dataset
            The ref dataset object used for setting ``self.ref_name_yrs``.
        time_selection : TimeSelection | None
            The optional time slice or season.
        """
        self.test_name_yrs = ds_test.get_name_yrs_attr(time_selection)
        self.ref_name_yrs = ds_ref.get_name_yrs_attr(time_selection)

    def _is_plevs_set(self):
        if (isinstance(self.plevs, np.ndarray) and not self.plevs.all()) or (
            not isinstance(self.plevs, np.ndarray) and not self.plevs
        ):
            return False

        return True

    def _run_diag(self) -> list[Any]:
        """Run the diagnostics for each set in the parameter.

        Additional CoreParameter (or CoreParameter sub-class) objects are derived
        from the CoreParameter `sets` attribute, hence this function returns a
        list of CoreParameter objects.

        This method loops over the parameter's diagnostic sets and attempts to
        import and call the related `run_diags()` function.

        Returns
        -------
        list[Any]
            The list of CoreParameter objects with results from the diagnostic run.
            NOTE: `Self` type is not yet supported by mypy.
        """
        results = []

        for set_name in self.sets:
            self.current_set = set_name
            # FIXME: This is a shortcut to importing `run_diag`, but can break
            # easily because the module driver is statically imported via string.
            # Instead, the import should be done more progammatically via
            # direct Python import.
            mod_str = "e3sm_diags.driver.{}_driver".format(set_name)

            # Check if there is a matching driver module for the `set_name`.
            try:
                module = importlib.import_module(mod_str)
            except ModuleNotFoundError as e:
                logger.error(f"'Error with set name {set_name}'", exc_info=e)
                continue

            # If the module exists, call the driver module's `run_diag` function.
            try:
                single_result = module.run_diag(self)
                results.append(single_result)
            except Exception:
                logger.exception(f"Error in {mod_str}", exc_info=True)

                if self.debug:
                    sys.exit()

        return results

    def _add_time_series_filepath_attr(
        self,
        data_type: Literal["test", "ref"],
        ds: xr.Dataset,
    ):
        """Add file path attributes to the parameter object.

        Parameters
        ----------
        data_type : Literal["test", "ref"]
            The type of data, either "test" or "ref".
        ds : xr.Dataset
            The dataset object containing the file path attribute.

        Raises
        ------
        ValueError
            If `data_type` is not "test" or "ref".
        """
        if data_type not in {"test", "ref"}:
            raise ValueError("data_type must be either 'test' or 'ref'.")

        filepath_attr = f"{data_type}_data_file_path"

        setattr(self, filepath_attr, getattr(ds, "file_path", "Unknown"))

    def _add_filepath_attr(
        self,
        data_type: Literal["test", "ref"],
        filepath: str | None = None,
    ):
        """Add file path attributes to the parameter object.

        Parameters
        ----------
        data_type : Literal["test", "ref"]
            The type of data, either "test" or "ref".
        filepath : str | None, optional
            The file path for climatology or time-slice data.

        Raises
        ------
        ValueError
            If `data_type` is not "test" or "ref".
        ValueError
            If `filepath` is not provided for climatology or time-slice data.
        """
        if data_type not in {"test", "ref"}:
            raise ValueError("data_type must be either 'test' or 'ref'.")

        filepath_attr = f"{data_type}_data_file_path"

        if not filepath:
            raise ValueError(
                "Filepath must be provided for the climatology or time-slice data."
            )

        setattr(self, filepath_attr, os.path.abspath(filepath))

    def _get_time_selection_to_use(
        self, require_one: bool = True
    ) -> tuple[Literal["time_slices", "seasons"], list[TimeSlice] | list[ClimoFreq]]:
        """
        Determine the time selection type and corresponding values.

        If ``time_slices`` are specified, they take precedence over ``seasons``.

        Parameters
        ----------
        require_one : bool, optional
            If True, ensures that at least one of `seasons` or `time_slices` is
            specified. Default is True.

        Returns
        -------
        tuple[Literal["time_slice", "seasons"], list[TimeSlice] | list[ClimoFreq]]
            A tuple containing the time selection type ("time_slices" or "seasons")
            and the corresponding list of values.

        Raises
        ------
        RuntimeError
            If neither `seasons` nor `time_slices` are specified when `require_one`
            is True.
        """
        if require_one and not (self.seasons or self.time_slices):
            raise RuntimeError(
                "Must specify either 'seasons' or 'time_slices'.\n"
                "- Use 'seasons' for climatological analysis (e.g., ['ANN', 'DJF']).\n"
                "- Use 'time_slices' for snapshot-based selection (e.g., ['0'], ['5'])."
            )

        # Time slices take precedence over seasons if both are specified.
        if self.time_slices:
            return "time_slices", self.time_slices

        return "seasons", self.seasons

    def __setattr__(self, name: str, value: Any) -> None:
        """Override setattr to ensure year attributes are padded when set."""
        if name in YEAR_ATTRIBUTES and value not in [None, ""]:
            # Validate and pad the year before setting the attribute
            value = pad_year(value)

        # Set the attribute using the superclass method
        super().__setattr__(name, value)
