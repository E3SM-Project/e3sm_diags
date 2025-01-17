"""This module stores the Dataset class, which is the primary class for I/O.

This Dataset class operates on `xr.Dataset` objects, which are created using
netCDF files. These `xr.Dataset` contain either the reference or test variable.
This variable can either be from a climatology file or a time series file.
If the variable is from a time series file, the climatology of the variable is
calculated. Reference and test variables can also be derived using other
variables from dataset files.
"""

from __future__ import annotations

import collections
import fnmatch
import glob
import os
import re
from datetime import datetime, timedelta
from typing import TYPE_CHECKING, Callable, Dict, List, Literal, Tuple

import dask
import pandas as pd
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.derivations import (
    DERIVED_VARIABLES,
    FUNC_NEEDS_TARGET_VAR,
    DerivedVariableMap,
    DerivedVariablesMap,
)
from e3sm_diags.driver import FRAC_REGION_KEYS, LAND_OCEAN_MASK_PATH
from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQS, ClimoFreq, climo
from e3sm_diags.driver.utils.regrid import HYBRID_SIGMA_KEYS
from e3sm_diags.logger import custom_logger

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter

logger = custom_logger(__name__)

# A constant variable that defines the pattern for time series filenames.
# Example: "ts_global_200001_200112.nc" (<VAR>_<SITE>_<TS_EXT_FILEPATTERN>)
TS_EXT_FILEPATTERN = r"_.{13}.nc"


# Additional variables to keep when subsetting.
HYBRID_VAR_KEYS = set(list(sum(HYBRID_SIGMA_KEYS.values(), ())))

# In some cases, lat and lon are stored as single point data variables rather
# than coordinates. These variables are kept when subsetting for downstream
# operations (e.g., arm_diags).
MISC_VARS = ["area", "areatotal2", "lat", "lon"]

# Seasons for model only data, used for matching filenames to construct
# filepaths.
MODEL_ONLY_SEASONS = [
    "ANN",
    "DJF",
    "MAM",
    "JJA",
    "SON",
    "01",
    "02",
    "03",
    "04",
    "05",
    "06",
    "07",
    "08",
    "09",
    "10",
    "11",
    "12",
]


def squeeze_time_dim(ds: xr.Dataset) -> xr.Dataset:
    """Squeeze single coordinate climatology time dimensions.

    For example, "ANN" averages over the year and collapses the time dim.
    Time bounds are also dropped if they exist.

    Parameters
    ----------
    ds : xr.Dataset
        The dataset with a time dimension

    Returns
    -------
    xr.Dataset
        The dataset with a time dimension.
    """
    try:
        time_dim = xc.get_dim_coords(ds, axis="T")
    except KeyError:
        time_dim = None

    if time_dim is not None and len(time_dim) == 1:
        ds = ds.squeeze(dim=time_dim.name)
        ds = ds.drop_vars(time_dim.name)

        bnds_key = time_dim.attrs.get("bounds")
        if bnds_key is not None and bnds_key in ds.data_vars.keys():
            ds = ds.drop_vars(bnds_key)

    return ds


class Dataset:
    def __init__(
        self,
        parameter: CoreParameter,
        data_type: Literal["ref", "test"],
    ):
        # The CoreParameter object with a list of parameters.
        self.parameter = parameter

        # The type of data for the Dataset object to store.
        self.data_type = data_type

        # The path, start year, and end year based on the dataset type.
        if self.data_type == "ref":
            self.root_path = self.parameter.reference_data_path
        elif self.data_type == "test":
            self.root_path = self.parameter.test_data_path
        else:
            raise ValueError(
                f"The `type` ({self.data_type}) for this Dataset object is invalid."
                "Valid options include 'ref' or 'test'."
            )

        # If the underlying data is a time series, set the `start_yr` and
        # `end_yr` attrs based on the data type (ref or test). Note, these attrs
        # are different for the `area_mean_time_series` parameter.
        if self.is_time_series:
            # FIXME: This conditional should not assume the first set is
            # area_mean_time_series. If area_mean_time_series is at another
            # index, this conditional is not False.
            if self.parameter.sets[0] in ["area_mean_time_series"]:
                self.start_yr = self.parameter.start_yr  # type: ignore
                self.end_yr = self.parameter.end_yr  # type: ignore
            elif self.data_type == "ref":
                self.start_yr = self.parameter.ref_start_yr  # type: ignore
                self.end_yr = self.parameter.ref_end_yr  # type: ignore
            elif self.data_type == "test":
                self.start_yr = self.parameter.test_start_yr  # type: ignore
                self.end_yr = self.parameter.test_end_yr  # type: ignore

        # The derived variables defined in E3SM Diags. If the `CoreParameter`
        # object contains additional user derived variables, they are added
        # to `self.derived_vars`.
        self.derived_vars_map = self._get_derived_vars_map()

        # Whether the data is sub-monthly or not.
        self.is_sub_monthly = False
        if self.parameter.sets[0] in ["diurnal_cycle", "arm_diags"]:
            self.is_sub_monthly = True

        # A boolean to keep track of whether the source variable(s) for a
        # target variable contains a wildcard ("?"). This is used to determine
        # whether to pass a list of wildcard variables as args to derived
        # variable function (True), or to unpack an expected number of variables
        # as args to a derived variable function (False).
        self.is_src_vars_wildcard = False

    @property
    def is_time_series(self):
        if (self.data_type == "ref" and self.parameter.ref_timeseries_input) or (
            self.data_type == "test" and self.parameter.test_timeseries_input
        ):
            return True
        else:
            return False

    @property
    def is_climo(self):
        return not self.is_time_series

    def _get_derived_vars_map(self) -> DerivedVariablesMap:
        """Get the defined derived variables.

        If the user-defined derived variables are in the input parameters,
        append parameters.derived_variables to the correct part of the derived
        variables dictionary.

        Returns
        -------
        DerivedVariablesMap
            A dictionary mapping the key of a derived variable to an ordered
            dictionary that maps a tuple of source variable(s) to a derivation
            function.
        """
        dvars: DerivedVariablesMap = DERIVED_VARIABLES.copy()
        user_dvars: DerivedVariablesMap = self.parameter.derived_variables

        # If the user-defined derived vars already exist, create a
        # new OrderedDict that combines the user-defined entries with the
        # existing ones in `e3sm_diags`. The user-defined entry should
        # be the highest priority and must be first in the OrderedDict.
        if user_dvars is not None:
            for key, ordered_dict in user_dvars.items():
                if key in dvars.keys():
                    dvars[key] = collections.OrderedDict(**ordered_dict, **dvars[key])  # type: ignore
                else:
                    dvars[key] = ordered_dict

        return dvars

    # Attribute related methods
    # --------------------------------------------------------------------------
    def get_name_yrs_attr(
        self,
        season: ClimoFreq | None = None,
        default_name: str | None = None,
    ) -> str:
        """Get the diagnostic name and 'yrs_averaged' attr as a single string.

        This method is used to update either `parameter.test_name_yrs` or
        `parameter.ref_name_yrs`, depending on `self.data_type`.

        If the dataset is contains a climatology, attempt to get "yrs_averaged"
        from the global attributes of the netCDF file. If this attribute cannot
        be retrieved, only return the diagnostic name.

        Parameters
        ----------
        season : CLIMO_FREQ | None, optional
            The climatology frequency, by default None.

        Returns
        -------
        str
            The name and years average string.
            Example: "historical_H1 (2000-2002)"

        Notes
        -----
        Replaces `e3sm_diags.driver.utils.general.get_name_and_yrs`
        """
        if self.data_type == "test":
            diag_name = self._get_test_name(default_name)
        elif self.data_type == "ref":
            diag_name = self._get_ref_name(default_name)

        if self.is_climo:
            if season is None:
                raise ValueError(
                    "A `season` argument must be supplied for climatology datasets "
                    "to try to get the global attribute 'yrs_averaged'."
                )

            yrs_averaged_attr = self._get_global_attr_from_climo_dataset(
                "yrs_averaged", season
            )

            if yrs_averaged_attr is None:
                return diag_name

        elif self.is_time_series:
            yrs_averaged_attr = f"{self.start_yr}-{self.end_yr}"

        return f"{diag_name} ({yrs_averaged_attr})"

    def _get_test_name(self, default_name: str | None = None) -> str:
        """Get the diagnostic test name.

        Returns
        -------
        str
            The diagnostic test name.

        Notes
        -----
        Replaces `e3sm_diags.driver.utils.general.get_name`
        """
        if self.parameter.short_test_name != "":
            return self.parameter.short_test_name
        elif self.parameter.test_name != "":
            return self.parameter.test_name
        else:
            if default_name is not None:
                return default_name

            # NOTE: This else statement is preserved from the previous CDAT
            # codebase to maintain the same behavior.
            if self.parameter.test_name == "":
                logger.warning("No test name was set, using empty string as test name.")

            return self.parameter.test_name

    def _get_ref_name(self, default_name: str | None = None) -> str:
        """Get the diagnostic reference name.

        Returns
        -------
        str
            The diagnostic reference name.

        Notes
        -----
        Replaces `e3sm_diags.driver.utils.general.get_name`
        """
        if self.parameter.short_ref_name != "":
            return self.parameter.short_ref_name
        elif self.parameter.reference_name != "":
            return self.parameter.reference_name
        else:
            # Covers cases such as streamflow which set the ref name to "GSIM".
            if default_name is not None:
                return default_name

            # NOTE: This else statement is preserved from the previous CDAT
            # codebase to maintain the same behavior.
            if self.parameter.ref_name == "":
                logger.warning(
                    "No reference name was set, using empty string as reference name."
                )

            return self.parameter.ref_name

    def _get_global_attr_from_climo_dataset(
        self, attr: str, season: ClimoFreq
    ) -> str | None:
        """Get the global attribute from the climo file based on the season.

        Parameters
        ----------
        attr : str
            The attribute to get (e.g., "Convention").
        season : CLIMO_FREQ
            The climatology frequency.

        Returns
        -------
        str | None
            The attribute string if it exists, otherwise None.
        """
        attr_val = None

        try:
            filepath = self._get_climo_filepath(season)
        except OSError:
            pass
        else:
            ds = self._open_climo_dataset(filepath)
            attr_val = ds.attrs.get(attr)

        return attr_val

    # --------------------------------------------------------------------------
    # Climatology related methods
    # --------------------------------------------------------------------------
    def get_climo_dataset(self, var: str, season: ClimoFreq) -> xr.Dataset:
        """Get the dataset containing the climatology variable.

        These variables can either be from the test data or reference data.
        If the variable is in a time series dataset, use the variable to
        calculate the climatology based on the selected frequency. If the
        variable is already a climatology variable, return the climatology
        dataset directly.

        Parameters
        ----------
        var : str
            The key of the climatology or time series variable to get the
            dataset for.
        season : CLIMO_FREQ, optional
            The season for the climatology.

        Returns
        -------
        xr.Dataset
            The dataset containing the climatology variable.

        Raises
        ------
        ValueError
            If the specified variable is not a valid string.
        ValueError
            If the specified season is not a valid string.
        """
        self.var = var

        if not isinstance(self.var, str) or self.var == "":
            raise ValueError("The `var` argument is not a valid string.")
        if not isinstance(season, str) or season not in CLIMO_FREQS:
            raise ValueError(
                "The `season` argument is not a valid string. Options include: "
                f"{CLIMO_FREQS}"
            )

        if self.is_time_series:
            ds = self.get_time_series_dataset(var)

            ds_climo = climo(ds, self.var, season).to_dataset()
            ds_climo = ds_climo.bounds.add_missing_bounds(axes=["X", "Y"])

            return ds_climo

        ds = self._get_climo_dataset(season)

        return ds

    def _get_climo_dataset(self, season: str) -> xr.Dataset:
        """Get the climatology dataset for the variable and season.

        Parameters
        ----------
        season : str
            The season for the climatology.

        Returns
        -------
        xr.Dataset
            The climatology dataset.

        Raises
        ------
        IOError
            If the variable was not found in the dataset or able to be derived
            using other datasets.
        """
        filepath = self._get_climo_filepath(season)
        logger.info(f"Opening climatology file: {filepath}")

        ds = self._open_climo_dataset(filepath)

        if self.var in self.derived_vars_map:
            ds = self._get_dataset_with_derived_climo_var(ds)
        elif self.var in ds.data_vars.keys():
            pass
        else:
            raise IOError(
                f"Variable '{self.var}' was not in the file '{filepath}', nor was "
                "it defined in the derived variables dictionary."
            )

        ds = squeeze_time_dim(ds)
        ds = self._subset_vars_and_load(ds, self.var)

        return ds

    def _open_climo_dataset(self, filepath: str) -> xr.Dataset:
        """Open a climatology dataset.

        Some climatology files have "time" as a scalar variable. If the scalar
        variable is a single integer instead of a 1D array with a length
        matching the equivalent dimension size, Xarray will `raise ValueError:
        dimension 'time' already exists as a scalar variable`. For this case,
        the "time" scalar variable is dropped when opening the dataset.

        If the scalar variable is dropped or climatology file only has a
        "time" dimension without coordinates, new "time" coordinates will be
        added to the dataset.

        Related issue: https://github.com/pydata/xarray/issues/1709

        Parameters
        ----------
        filepath : str
            The path to a climatology dataset.

        Returns
        -------
        xr.Dataset
            The climatology dataset.

        Raises
        ------
        ValueError
            Raised for all ValueErrors other than "dimension 'time' already
            exists as a scalar variable".
        """
        # Time coordinates are decoded because there might be cases where
        # a multi-file climatology dataset has different units between files
        # but raw encoded time values overlap. Decoding with Xarray allows
        # concatenation of datasets with this issue (e.g., `area_cycle_zonal_mean`
        # set with the MERRA2_Aerosols climatology datasets).
        # NOTE: This GitHub issue explains why the "coords" and "compat" args
        # are defined as they are below: https://github.com/xCDAT/xcdat/issues/641
        args = {
            "paths": filepath,
            "decode_times": True,
            "add_bounds": ["X", "Y", "Z"],
            "coords": "minimal",
            "compat": "override",
        }

        try:
            ds = xc.open_mfdataset(**args)
        except ValueError as e:  # pragma: no cover
            # FIXME: Need to fix the test that covers this code block.
            msg = str(e)

            if "dimension 'time' already exists as a scalar variable" in msg:
                ds = xc.open_mfdataset(**args, drop_variables=["time"])
            else:
                raise ValueError(msg) from e

        if "time" not in ds.coords:
            ds["time"] = xr.DataArray(
                name="time",
                dims=["time"],
                data=[0],
                attrs={"axis": "T", "standard_name": "time"},
            )

        return ds

    def _get_climo_filepath(self, season: str) -> str:
        """Return the path to the climatology file.

        There are three patterns for matching a file, with the first match
        being returned if any match is found:

        1. Using the reference/test file parameters if they are set (`ref_file`,
           `test_file`).
           - {reference_data_path}/{ref_file}
           - {test_data_path}/{test_file}
        2. Using the reference/test name and season.
           - {reference_data_path}/{ref_name}_{season}.nc
           - {test_data_path}/{test_name}_{season}.nc
        3. Using the reference or test name as a nested directory with the same
           name as the filename with a season.
           - General match pattern:
             - {reference_data_path}/{ref_name}/{ref_name}_{season}.nc
             - {test_data_path}/{test_name}/{test_name}_{season}.nc
           - Patern for model-only data for season in "ANN" "DJF", "MAM", "JJA",
             or "SON":
             - {reference_data_path}/{ref_name}/{ref_name}.*{season}.*.nc
             - {test_data_path}/{test_name}/{test_name}.*{season}.*.nc

        Parameters
        ----------
        season : str
            The season for the climatology.

        Returns
        -------
        str
            The path to the climatology file.
        """
        # First pattern attempt.
        filepath = self._get_climo_filepath_with_params()

        # Second and third pattern attempts.
        if filepath is None:
            if self.data_type == "ref":
                filename = self.parameter.ref_name
            elif self.data_type == "test":
                filename = self.parameter.test_name

            if season == "ANNUALCYCLE":
                filepath = self._find_climo_filepath(filename, "01")

                # find the path for 12 monthly mean files
                if filepath is not None:
                    filename_01 = filepath.split("/")[-1]
                    filepath = filepath.replace(
                        # f"{filename_01}", f"{filename}_[0-1][0-9]_*_*climo.nc"
                        # AOD_550 dataset has pattern AOD_550_01_climo.nc, other dataset has e.g ERA5_ANN_197901_201912_climo.nc
                        f"{filename_01}",
                        f"{filename}_[0-1][0-9]_*climo.nc",
                    )
            else:
                filepath = self._find_climo_filepath(filename, season)

            # If absolutely no filename was found, then raise an error.
            if filepath is None:
                raise IOError(
                    f"No file found for '{filename}' and '{season}' in {self.root_path}"
                )

        return filepath

    def _get_climo_filepath_with_params(self) -> str | None:
        """Get the climatology filepath using parameters.

        Returns
        -------
        str | None
            The filepath using the `ref_file` or `test_file`  parameter if they
            are set.
        """
        filepath = None

        if self.data_type == "ref":
            if self.parameter.ref_file != "":
                filepath = os.path.join(self.root_path, self.parameter.ref_file)

        elif self.data_type == "test":
            if hasattr(self.parameter, "test_file"):
                filepath = os.path.join(self.root_path, self.parameter.test_file)

        return filepath

    def _find_climo_filepath(self, filename: str, season: str) -> str | None:
        """Find the climatology filepath for the variable.

        Parameters
        ----------
        filename : str
            The filename for the climatology variable.
        season : str
            The season for climatology.

        Returns
        -------
        str | None
            The filepath for the climatology variable.
        """
        # First attempt: try to find the climatology file based on season.
        # Example: {path}/{filename}_{season}.nc
        filepath = self._find_climo_filepath_with_season(
            self.root_path, filename, season
        )

        # Second attempt: try looking for the file nested in a folder, based on
        # the test_name.
        # Example: {path}/{filename}/{filename}_{season}.nc
        #           data_path/some_file/some_file_ANN.nc
        if filepath is None:
            nested_root_path = os.path.join(self.root_path, filename)

            if os.path.exists(nested_root_path):
                filepath = self._find_climo_filepath_with_season(
                    nested_root_path, filename, season
                )

        return filepath

    def _find_climo_filepath_with_season(
        self, root_path: str, filename: str, season: str
    ) -> str | None:
        """Find climatology filepath with a root path, filename, and season.

        Parameters
        ----------
        root_path : str
            The root path containing `.nc` files. The `.nc` files can be nested
            in sub-directories within the root path.
        filename : str
            The filename for the climatology variable.
        season : str
            The season for climatology.

        Returns
        -------
        str | None
            The climatology filepath based on season, if it exists.
        """
        files_in_dir = sorted(os.listdir(root_path))

        # If the filename is followed by _<SEASON>.
        for file in files_in_dir:
            if file.startswith(filename + "_" + season):
                return os.path.join(root_path, file)

        # For model only data, the <SEASON> string can by anywhere in the
        # filename. This is a more general pattern for model only data.
        if season in MODEL_ONLY_SEASONS:
            for file in files_in_dir:
                if file.startswith(filename) and season in file:
                    return os.path.join(root_path, file)

        return None

    def _get_dataset_with_derived_climo_var(self, ds: xr.Dataset) -> xr.Dataset:
        """Get the dataset containing the derived variable (`self.var`).

        Parameters
        ----------
        ds: xr.Dataset
            The climatology dataset, which should contain the source variables
            for deriving the target variable.

        Returns
        -------
        xr.Dataset
            The dataset with the derived variable.

        Raises
        ------
        IOError
            If the datasets for the target variable and source variables were
            not found in the data directory, or the target variable cannot be
            found directly in the dataset.
        """
        # An OrderedDict mapping possible source variables to the function
        # for deriving the variable of interest.
        # Example: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}
        target_var = self.var
        target_var_map = self.derived_vars_map[target_var]

        # Get the first valid source variables and its derivation function.
        # The source variables are checked to exist in the dataset object
        # and the derivation function is used to derive the target variable.
        # Example:
        #   For target variable "PRECT": {('PRECC', 'PRECL'): func}
        matching_target_var_map = self._get_matching_climo_src_vars(ds, target_var_map)

        if matching_target_var_map is not None:
            # NOTE: Since there's only one set of vars, we get the first and only set
            # of vars from the derived variable dictionary.
            # 1. Get the derivation function.
            derivation_func = list(matching_target_var_map.values())[0]

            # 2. Get the derivation function arguments using source variable keys.
            # Example: [xr.DataArray(name="PRECC",...), xr.DataArray(name="PRECL",...)]
            src_var_keys = list(matching_target_var_map.keys())[0]

            logger.info(
                f"Deriving the climatology variable using the source variables: {src_var_keys}"
            )
            ds_sub = squeeze_time_dim(ds)
            ds_sub = self._subset_vars_and_load(ds_sub, list(src_var_keys))

            # 3. Use the derivation function to derive the variable.
            ds_derived = self._get_dataset_with_derivation_func(
                ds_sub, derivation_func, src_var_keys, target_var
            )

            return ds_derived

        # None of the entries in the derived variables dictionary worked,
        # so try to get the variable directly from the dataset.
        if target_var in ds.data_vars.keys():
            return ds

        raise IOError(
            f"The dataset file has no matching source variables for {target_var} and "
            f"could not be derived using {list(target_var_map.keys())}."
        )

    def _get_matching_climo_src_vars(
        self,
        dataset: xr.Dataset,
        target_variable_map: DerivedVariableMap,
    ) -> DerivedVariableMap | None:
        """Get the matching climatology source vars based on the target variable.

        Parameters
        ----------
        dataset : xr.Dataset
            The dataset containing the source variables.
        target_var : str
            The target variable to derive.
        target_var_map : TARGET_VARIABLE_MAP
            An ordered dictionary mapping the target variable's source variables
            to their derivation functions.

        Returns
        -------
        DerivedVariableMap | None
            The optional matching dictionary with the key being the source
            variables and the value being the derivation function.
        """
        self.is_src_vars_wildcard = False

        vars_in_file = set(dataset.data_vars.keys())

        # Example: [('pr',), ('PRECC', 'PRECL')]
        possible_vars = list(target_variable_map.keys())

        # Try to get the var using entries from the dictionary.
        for var_tuple in possible_vars:
            var_list = list(var_tuple).copy()

            for vars in var_tuple:
                # Add support for wild card `?` in variable strings
                # Example: ('bc_a?DDF', 'bc_c?DDF')
                if "?" in vars:
                    var_list += fnmatch.filter(list(vars_in_file), vars)
                    var_list.remove(vars)
                    self.is_src_vars_wildcard = True

            if vars_in_file.issuperset(tuple(var_list)):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {tuple(var_list): target_variable_map[var_tuple]}

        return None

    # --------------------------------------------------------------------------
    # Time series related methods
    # --------------------------------------------------------------------------
    def get_time_series_dataset(
        self, var: str, single_point: bool = False
    ) -> xr.Dataset:
        """Get variables from time series datasets.

        Variables must exist in the time series files. These variables can
        either be from the test data or reference data.

        Parameters
        ----------
        var : str
            The key of the time series variable to get the dataset for.
        single_point : bool, optional
            Single point indicating the data is sub monthly, by default False.
            If False, center the time coordinates using time bounds.

        Returns
        -------
        xr.Dataset
            The time series Dataset.

        Raises
        ------
        ValueError
            If the dataset is not a time series.
        ValueError
            If the `var` argument is not a string or an empty string.
        IOError
            If the variable does not have a file in the specified directory
            and it was not defined in the derived variables dictionary.
        """
        self.var = var

        if not self.is_time_series:
            raise ValueError("You can only use this function with time series data.")

        if not isinstance(self.var, str) or self.var == "":
            raise ValueError("The `var` argument is not a valid string.")

        if self.var in self.derived_vars_map:
            ds = self._get_dataset_with_derived_ts_var()
        else:
            ds = self._get_time_series_dataset_obj(self.var)

        if single_point:
            self.is_sub_monthly = True

        if not self.is_sub_monthly:
            ds = self._center_time_for_non_submonthly_data(ds)

        return ds

    def _get_dataset_with_derived_ts_var(self) -> xr.Dataset:
        """Get the dataset containing the derived time series variable.

        Returns
        -------
        xr.Dataset
            The dataset with the derived time series variable.
        """
        # An OrderedDict mapping possible source variables to the function
        # for deriving the variable of interest.
        # Example: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}
        target_var = self.var
        target_var_map = self.derived_vars_map[target_var]

        # Get the first valid source variables and its derivation function.
        # The source variables are checked to exist in the dataset object
        # and the derivation function is used to derive the target variable.
        # Example:
        #   For target variable "PRECT": {('PRECC', 'PRECL'): func}
        matching_target_var_map = self._get_matching_time_series_src_vars(
            self.root_path, target_var_map
        )

        # NOTE: Since there's only one set of vars, we get the first and only set
        # of vars from the derived variable dictionary.
        # 1. Get the derivation function.
        derivation_func = list(matching_target_var_map.values())[0]

        # 2. Get the derivation function arguments using source variable keys.
        # Example: [xr.DataArray(name="PRECC",...), xr.DataArray(name="PRECL",...)]
        # Unlike the climatology dataset, the source variables for time series
        # data can be found in multiple datasets so a single xr.Dataset object
        # is returned containing all of them.
        src_var_keys = list(matching_target_var_map.keys())[0]

        logger.info(
            f"Deriving the time series variable using the source variables: {src_var_keys}"
        )

        ds = self._get_dataset_with_source_vars(src_var_keys)

        # 3. Use the derivation function to derive the variable.
        ds_final = self._get_dataset_with_derivation_func(
            ds, derivation_func, src_var_keys, target_var
        )

        return ds_final

    def _get_dataset_with_derivation_func(
        self,
        ds: xr.Dataset,
        func: Callable,
        src_var_keys: Tuple[str, ...],
        target_var_key: str,
    ) -> xr.Dataset:
        """
        Get the dataset with the target variable using the derivation function
        and source variables.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset with source variables used for deriving the target
            variable.
        func : Callable
            The derivation function that uses the source variables to derive
            the target variables.
        src_var_keys : Tuple[str, ...]
            The source variable keys.
        target_var_key : str
            The target variable key.

        Returns
        -------
        xr.Dataset
            The dataset with the derived target variable.
        """
        func_args = [ds[var].copy() for var in src_var_keys]

        if func in FUNC_NEEDS_TARGET_VAR:
            func_args = [target_var_key] + func_args  # type: ignore # pragma: nocover

        # If the target variable key contains a wildcard, there are can be
        # N number of function arguments. For the cases, we need to pass the
        # entire list of DataArrays to the derived function in order to
        # properly maintain attributes (e.g., "units"). For example,
        # "bc_a?_CLXF" uses the Python `sum()` function, which sums all of the
        # DataArrays by iterating over them and add all elements by index (rather
        # than a typical sum reduction across an axis). Using `.sum()` without
        # using `with xr.set_options(keep_attrs=True)` will result in attributes
        # being dropped, resulting in potential downstream issues such as unit
        # conversions which require the "units" attribute.
        if self.is_src_vars_wildcard:
            derived_var = func(func_args)  # pragma: nocover
        else:
            derived_var = func(*func_args)  # pragma: nocover

        # Derive the target variable, then align the dataset dimensions to it.
        # Dimensional alignment is necessary for cases where the target variable
        # has been subsetted (e.g., `cosp_histogram_standardize()`). `xr.align`
        # returns both objects and element 0 is the xr.Dataset that is needed.
        ds_final = xr.align(ds.copy(), derived_var)[0]
        ds_final[target_var_key] = derived_var

        return ds_final

    def _get_matching_time_series_src_vars(
        self, path: str, target_var_map: DerivedVariableMap
    ) -> Dict[Tuple[str, ...], Callable]:
        """Get the matching time series source vars based on the target variable.

        Parameters
        ----------
        path: str
            The path containing the dataset(s).
        target_var_map : DerivedVariableMap
            An ordered dictionary for a target variable that maps a tuple of
            source variable(s) to a derivation function.

        Returns
        -------
        DerivedVariableMap
            The matching dictionary with the key being the source variable(s)
            and the value being the derivation function.

        Raises
        ------
        IOError
            If the datasets for the target variable and source variables were
            not found in the data directory.
        """
        # Example: [('pr',), ('PRECC', 'PRECL')]
        possible_vars = list(target_var_map.keys())

        # Loop over the tuples of possible source variable and try to get
        # the matching derived variables dictionary if the files exist in the
        # time series filepath.
        for tuple_of_vars in possible_vars:
            all_vars_found = all(
                self._get_time_series_filepaths(path, var) is not None
                for var in tuple_of_vars
            )

            if all_vars_found:
                return {tuple_of_vars: target_var_map[tuple_of_vars]}

        # None of the entries in the derived variables dictionary are valid,
        # so try to get the dataset for the variable directly.
        # Example file name: {var}_{start_yr}01_{end_yr}12.nc.
        if self._get_time_series_filepaths(path, self.var) is not None:
            return {(self.var,): lambda x: x}

        raise IOError(
            f"No files found for target variable {self.var} or derived variables "
            f"({possible_vars}) in {path}."
        )

    def _get_dataset_with_source_vars(self, vars_to_get: Tuple[str, ...]) -> xr.Dataset:
        """Get the variables from datasets in the specified path.

        Parameters
        ----------
        path : str
            The path to the datasets.
        vars_to_get: Tuple[str]
            The source variables used to derive the target variable.

        Returns
        -------
        xr.Dataset
            The dataset with the source variables.
        """
        datasets = []

        for var in vars_to_get:
            ds = self._get_time_series_dataset_obj(var)
            datasets.append(ds)

        ds = xr.merge(datasets)
        ds = squeeze_time_dim(ds)

        return ds

    def _get_time_series_dataset_obj(self, var) -> xr.Dataset:
        """Get the time series dataset for a variable.

        This method also parses the start and end time from the dataset filename
        to subset the dataset.

        Returns
        -------
        xr.Dataset
            The dataset for the variable.
        """
        filepaths = self._get_time_series_filepaths(self.root_path, var)
        logger.info(f"Opening time series files: {filepaths}")

        if filepaths is None:
            raise IOError(
                f"No time series `.nc` file was found for '{var}' in '{self.root_path}'"
            )

        ds = xc.open_mfdataset(
            filepaths,
            add_bounds=["X", "Y", "T"],
            decode_times=True,
            use_cftime=True,
            coords="minimal",
            compat="override",
        )
        ds_subset = self._subset_time_series_dataset(ds, var)

        return ds_subset

    def _get_time_series_filepaths(
        self, root_path: str, var_key: str
    ) -> List[str] | None:
        """Get the matching variable time series filepaths.

        This method globs the specified path for all `*.nc` files and attempts
        to find the matching time series filepath(s) for the specified variable.

        Example matching filenames.
            - {var}_{start_yr}01_{end_yr}12.nc
            - {self.parameters.ref_name}/{var}_{start_yr}01_{end_yr}12.nc

        Parameters
        ----------
        root_path : str
            The root path containing `.nc` files. The `.nc` files can be nested
            in sub-directories within the root path.
        var_key : str
            The variable key used to find the time series file.

        Returns
        -------
        List[str]
            A list of matching filepaths for the variable. If no match is found,
            None is returned.
        """
        # The filename pattern for matching using regex.
        if self.parameter.sets[0] in ["arm_diags"]:
            # Example: "ts_global_200001_200112.nc"
            site = getattr(self.parameter, "regions", "")
            filename_pattern = var_key + "_" + site[0] + TS_EXT_FILEPATTERN
        else:
            # Example: "ts_200001_200112.nc"
            filename_pattern = var_key + TS_EXT_FILEPATTERN

        # First pattern example: {path}/ts_200001_200112.nc"
        matches = self._get_matches(root_path, filename_pattern)

        # If no matches were found with the first pattern, try the second
        # pattern using ref_name.
        # Second pattern example: {path}/{ref_name}/ts_200001_200112.nc"
        ref_name = getattr(self.parameter, "ref_name", None)
        if len(matches) == 0 and ref_name is not None:
            matches = self._get_matches(root_path, filename_pattern, ref_name)

        if len(matches) == 0:
            return None

        return matches

    def _get_matches(
        self, root_path: str, filename_pattern: str, ref_name: str | None = None
    ) -> List[str]:
        """Get the matching filepaths based on the glob path and pattern.

        Parameters
        ----------
        root_path : str
            The root path to search for files.
        filepath_pattern : str
            The regex pattern to match filepaths.
            For example, "RIVER_DISCHARGE_OVER_LAND_LIQ_.{13}.nc".
        ref_name : str | None
            The directory name storing references files, by default None.

        Returns
        -------
        List[str]
            A list of matching filepaths.
        """
        if ref_name is None:
            glob_dir = root_path
            filepath_pattern = os.path.join(root_path, filename_pattern)
        else:
            glob_dir = os.path.join(root_path, ref_name)
            filepath_pattern = os.path.join(root_path, ref_name, filename_pattern)

        glob_path = os.path.join(glob_dir, "**", "*.nc")
        filepaths = glob.glob(glob_path, recursive=True)
        filepaths = sorted(filepaths)
        matches = [f for f in filepaths if re.search(filepath_pattern, f)]

        return matches

    def _subset_time_series_dataset(self, ds: xr.Dataset, var: str) -> xr.Dataset:
        """Subset the time series dataset.

        This method subsets the variables in the dataset and loads the data
        into memory, then subsets on the time slice based on the specified
        files.

        Parameters
        ----------
        ds : xr.Dataset
            The time series dataset.
        var : str
            The main variable to keep.

        Returns
        -------
        xr.Dataset
            The subsetted time series dataset.

        """
        time_slice = self._get_time_slice(ds)
        ds_sub = ds.sel(time=time_slice).squeeze()

        if self.is_sub_monthly:
            ds_sub = self._exclude_sub_monthly_coord_spanning_year(ds_sub)

        ds_sub = self._subset_vars_and_load(ds_sub, var)

        return ds_sub

    def _get_time_slice(self, ds: xr.Dataset) -> slice:
        """Get time slice to subset a dataset.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset.
        Returns
        -------
        slice
            A slice object with a start and end time in the format "YYYY-MM-DD".

        Raises
        ------
        ValueError
            If invalid date range specified for test/reference time series data.
        """
        start_yr_int, end_yr_int = int(self.start_yr), int(self.end_yr)
        var_start_year, var_end_year = self._extract_var_start_and_end_years(ds)

        if start_yr_int < var_start_year:
            raise ValueError(
                "Invalid year range specified for test/reference time series data: "
                f"start_year ({start_yr_int}) < var_start_yr ({var_start_year})."
            )
        elif end_yr_int > var_end_year:
            raise ValueError(
                "Invalid year range specified for test/reference time series data: "
                f"end_year ({end_yr_int}) > var_end_yr ({var_end_year})."
            )

        start_yr_str = str(start_yr_int).zfill(4)
        end_yr_str = str(end_yr_int).zfill(4)

        if self.is_sub_monthly:
            start_time = f"{start_yr_str}-01-01"

            end_yr_str = str(int(end_yr_str) + 1).zfill(4)
            end_time = f"{end_yr_str}-01-01"
        else:
            start_time = self._get_slice_with_bounds(ds, start_yr_str, "start")
            end_time = self._get_slice_with_bounds(ds, end_yr_str, "end")

        return slice(start_time, end_time)

    def _extract_var_start_and_end_years(self, ds: xr.Dataset) -> Tuple[int, int]:
        """Extract the start and end years from the time coordinates.

        If the last time coordinate starts in January, subtract one year from
        the end year to get the correct end year which should align with the
        end year from the filepaths.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset with time coordinates.

        Returns
        -------
        Tuple[int, int]
            The start and end years.
        """
        time_coords = xc.get_dim_coords(ds, axis="T")
        start_year = time_coords[0].dt.year
        end_year = time_coords[-1].dt.year

        if time_coords[-1].dt.month == 1:
            end_year -= 1

        return start_year.item(), end_year.item()

    def _get_slice_with_bounds(
        self, ds: xr.Dataset, year_str: str, slice_type: Literal["start", "end"]
    ) -> str:
        """Get the time slice for non-submonthly data using bounds if needed.

        For example, let's say we have non-submonthly time coordinates with a
        start time and end time of ("2011-01-01", "2014-01-01"). We pre-define
        a time slice of ("2011-01-15", "2013-12-15"). However slicing with
        an end time of "2013-12-15" will exclude the last coordinate
        "2014-01-01". To rectify this situation, we use the time delta between
        time bounds to extend the end time slice to "2014-01-15".

          1. Get the time delta between bound values ["2013-12-15",
             "2014-01-15"], which is one month.
          2. Add the time delta to the old end time of "2013-12-15", which
             results in a new end time of "2014-01-15"
             - Time delta is subtracted if start time slice needs to be
             adjusted.
          3. Now slice the time coordinates using ("2011-01-15", "2014-01-15").
             Xarray will now correctly correctly subset to include the last
             coordinate value of "2014-01-01" using this time slice.
        Parameters
        ----------
        ds : xr.Dataset
            The dataset.
        year_str : str
            The year string for the slice.
        slice_type : Literal["start", "end"]
            The slice type, start or end.

        Returns
        -------
        str
            The new time slice type.

        Notes
        -----
        This function replicates the cdms2/cdutil "ccb" slice flag used for
        subsetting. "ccb" only allows the right side to be closed. It will get
        the difference between bounds values and adjust the subsetted start
        and end time coordinates depending on whether they fall within the
        slice or not.
        """
        time_bounds = ds.bounds.get_bounds(axis="T")
        time_delta = self._get_time_bounds_delta(time_bounds)
        time_coords = xc.get_dim_coords(ds, axis="T")
        actual_day, actual_month = (
            time_coords[0].dt.day.item(),
            time_coords[0].dt.month.item(),
        )

        if slice_type == "start":
            stop = f"{year_str}-01-15"

            if actual_day < 15 and actual_month == 1:
                stop = self._adjust_slice_str(stop, time_delta, add=False)
        else:
            stop = f"{year_str}-12-15"

            if actual_day > 15 or actual_month > 1:
                stop = self._adjust_slice_str(stop, time_delta, add=True)

        return stop

    def _adjust_slice_str(self, slice_str: str, delta: timedelta, add: bool) -> str:
        """Adjusts a date string by a given time delta.

        Parameters
        ----------
        slice_str : str
            The date string to be adjusted, in the format "%Y-%m-%d".
        delta : timedelta
            The time delta by which to adjust the date.
        add : bool
            If True, the delta is added to the date; if False, the delta is
            subtracted.

        Returns
        -------
        str
            The adjusted date string in ISO format.
        """
        slice_dt = datetime.strptime(slice_str, "%Y-%m-%d")
        slice_dt_new = slice_dt + delta if add else slice_dt - delta

        slice_str_new = self._convert_new_stop_pt_to_iso_format(slice_dt_new)

        return slice_str_new

    def _get_time_bounds_delta(self, time_bnds: xr.DataArray) -> timedelta:
        """Get the time delta between bounds values.

        Parameters
        ----------
        time_bnds : xr.DataArray
            The time bounds.

        Returns
        -------
        timedelta
            The time delta.
        """
        time_delta = time_bnds[0][-1] - time_bnds[0][0]

        # The source dataset object has not been loaded into memory yet, so
        # load the time delta into memory using the sync scheduler to avoid
        # hanging from conflicting schedulers.
        if isinstance(time_delta.data, dask.array.core.Array):
            time_delta.load(scheduler="sync")

        time_delta_py = pd.to_timedelta(time_delta.values).to_pytimedelta()

        return time_delta_py

    def _convert_new_stop_pt_to_iso_format(self, new_stop: datetime) -> str:
        """Convert the new stop point ISO-8061 formatted string.

        For example, "2012-12-15" and "0051-12-01". Otherwise, Xarray will
        raise `ValueError: no ISO-8601 or cftime-string-like match for string:`

        Parameters
        ----------
        new_stop : datetime
            The new stop point.

        Returns
        -------
        str
            The new stop point as an ISO-8061 formatted string.
        """
        year = str(new_stop.year).zfill(4)
        month = str(new_stop.month).zfill(2)
        day = str(new_stop.day).zfill(2)

        return f"{year}-{month}-{day}"

    def _exclude_sub_monthly_coord_spanning_year(
        self, ds_subset: xr.Dataset
    ) -> xr.Dataset:
        """
        Exclude the last time coordinate for sub-monthly data if it extends into
        the next year.

        Excluding end time coordinates that extend to the next year is
        necessary because downstream operations such as annual cycle climatology
        should consist of data for full years for accurate calculations.

        For example, if the time slice is ("0001-01-01", "0002-01-01") and
        the last time coordinate is:
            * "0001-12-31" -> don't exclude
            * "0002-01-01" -> exclude

        Parameters
        ----------
        ds_subset : xr.Dataset
            The subsetted dataset.

        Returns
        -------
        xr.Dataset
            The dataset with the last time coordinate excluded if necessary.

        Notes
        -----
        This function replicates the CDAT cdms2 "co" slice flag (close, open).
        """
        time_dim = xc.get_dim_keys(ds_subset, axis="T")
        time_values = ds_subset[time_dim]
        last_time_year = time_values[-1].dt.year.item()
        second_last_time_year = time_values[-2].dt.year.item()

        if last_time_year > second_last_time_year:
            ds_subset = ds_subset.isel(time=slice(0, -1))

        return ds_subset

    def _center_time_for_non_submonthly_data(self, ds: xr.Dataset) -> xr.Dataset:
        """Center time coordinates using bounds for non-submonthly data.

        This is important for data where the absolute time doesn't fall in the
        middle of the time interval, such as E3SM native format data where
        the time was recorded at the end of each time bounds.

        Parameters
        ----------
        ds : xr.Dataset
            The dataset.

        Returns
        -------
        ds : xr.Dataset
            The dataset with centered time coordinates.
        """
        try:
            time_dim = xc.get_dim_keys(ds, axis="T")
        except (ValueError, KeyError):
            time_dim = None

        if time_dim is not None:
            return xc.center_times(ds)

        return ds

    def _get_land_sea_mask(self, season: str) -> xr.Dataset:
        """Get the land sea mask from the dataset or use the default file.

        Land sea mask variables are time invariant which means the time
        dimension will be squeezed and dropped from the final xr.Dataset
        output since it is not needed.

        Parameters
        ----------
        season : str
            The season to subset on.

        Returns
        -------
        xr.Dataset
            The xr.Dataset object containing the land sea mask variables
            "LANDFRAC" and "OCNFRAC".
        """
        ds_mask = self._get_land_sea_mask_dataset(season)

        if ds_mask is None:
            logger.info("No land sea mask datasets were found for the given season.")
            ds_mask = self._get_default_land_sea_mask_dataset()

        return ds_mask

    def _get_land_sea_mask_dataset(self, season: str) -> xr.Dataset | None:
        """Get the land sea mask dataset for the given season.

        Parameters
        ----------
        season : str
            The season to subset on.

        Returns
        -------
        xr.Dataset | None
            The land sea mask dataset for the given season, or None if not
            found.
        """
        land_keys = FRAC_REGION_KEYS["land"]
        ocn_keys = FRAC_REGION_KEYS["ocean"]

        datasets = []

        for land_key, ocn_key in zip(land_keys, ocn_keys, strict=False):
            try:
                ds_land = self.get_climo_dataset(land_key, season)  # type: ignore
                ds_ocn = self.get_climo_dataset(ocn_key, season)  # type: ignore
            except IOError:
                pass
            else:
                datasets.append(ds_land)
                datasets.append(ds_ocn)

        if len(datasets) == 2:
            return xr.merge(datasets)

        return None

    def _get_default_land_sea_mask_dataset(self) -> xr.Dataset:
        """Get the default land sea mask dataset.

        Returns
        -------
        xr.Dataset
            The default land sea mask dataset.
        """
        logger.info(f"Using default land sea mask located at `{LAND_OCEAN_MASK_PATH}`.")

        ds_mask = xr.open_dataset(LAND_OCEAN_MASK_PATH)
        ds_mask = squeeze_time_dim(ds_mask)

        return ds_mask

    def _subset_vars_and_load(self, ds: xr.Dataset, var: str | List[str]) -> xr.Dataset:
        """Subset for variables needed for processing and load into memory.

        Subsetting the dataset reduces its memory footprint. Loading is
        necessary because there seems to be an issue with `open_mfdataset()`
        and using the multiprocessing scheduler defined in e3sm_diags,
        resulting in timeouts and resource locking. To avoid this, we load the
        multi-file dataset into memory before performing downstream operations.

        Source: https://github.com/pydata/xarray/issues/3781

        Parameters
        ----------
        ds : xr.Dataset
            The dataset.
        var : str | List[str]
            The variable or variables to subset on.

        Returns
        -------
        xr.Dataset
            The dataset subsetted and loaded into memory.
        """
        # slat and slon are lat lon pair for staggered FV grid included in
        # remapped files.
        if "slat" in ds.dims:
            ds = ds.drop_dims(["slat", "slon"])

        all_vars_keys = list(ds.data_vars.keys())
        keep_vars = [
            var
            for var in all_vars_keys
            if "bnd" in var
            or "bounds" in var
            or var in HYBRID_VAR_KEYS
            or var in MISC_VARS
        ]

        if isinstance(var, str):
            var = [var]

        ds = ds[var + keep_vars]

        ds.load(scheduler="sync")

        return ds
