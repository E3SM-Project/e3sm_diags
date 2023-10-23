"""This module stores the Dataset class, which is the primary class for I/O.

NOTE: Replaces `e3sm_diags.driver.utils.dataset`.

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
from typing import TYPE_CHECKING, Callable, Dict, Literal, Tuple

import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.derivations import (
    DERIVED_VARIABLES,
    DerivedVariableMap,
    DerivedVariablesMap,
)
from e3sm_diags.driver import LAND_FRAC_KEY, LAND_OCEAN_MASK_PATH, OCEAN_FRAC_KEY
from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQ, CLIMO_FREQS, climo
from e3sm_diags.logger import custom_logger

if TYPE_CHECKING:
    from e3sm_diags.parameter.core_parameter import CoreParameter


logger = custom_logger(__name__)

# A constant variable that defines the pattern for time series filenames.
# Example: "ts_global_200001_200112.nc" (<VAR>_<SITE>_<TS_EXT_FILEPATTERN>)
TS_EXT_FILEPATTERN = r"_.{13}.nc"


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

    @property
    def is_time_series(self):
        if self.parameter.ref_timeseries_input or self.parameter.test_timeseries_input:
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
        user_dvars: DerivedVariablesMap = getattr(self.parameter, "derived_variables")

        # If the user-defined derived vars already exist, create a
        # new OrderedDict that combines the user-defined entries with the
        # existing ones in `e3sm_diags`. The user-defined entry should
        # be the highest priority and must be first in the OrderedDict.
        if user_dvars is not None:
            for key, ordered_dict in user_dvars.items():
                if key in dvars.keys():
                    dvars[key] = collections.OrderedDict(**ordered_dict, **dvars[key])
                else:
                    dvars[key] = ordered_dict

        return dvars

    # Attribute related methods
    # --------------------------------------------------------------------------
    def get_name_yrs_attr(self, season: CLIMO_FREQ | None = None) -> str:
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
            diag_name = self._get_test_name()
        elif self.data_type == "ref":
            diag_name = self._get_ref_name()

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

    def _get_test_name(self) -> str:
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

        raise AttributeError(
            "Either `parameter.short_test_name` or `parameter.test_name attributes` "
            "must be set to get the name and years attribute for test datasets."
        )

    def _get_ref_name(self) -> str:
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
        elif self.parameter.ref_name != "":
            return self.parameter.ref_name

        raise AttributeError(
            "Either `parameter.short_ref_name`, `parameter.reference_name`, or "
            "`parameter.ref_name` must be set to get the name and years attribute for "
            "reference datasets."
        )

        return self.parameter.ref_name

    def _get_global_attr_from_climo_dataset(
        self, attr: str, season: CLIMO_FREQ
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
        filepath = self._get_climo_filepath(season)

        ds = xr.open_dataset(filepath)
        attr_val = ds.attrs.get(attr)

        return attr_val

    # --------------------------------------------------------------------------
    # Climatology related methods
    # --------------------------------------------------------------------------
    def get_ref_climo_dataset(
        self, var_key: str, season: CLIMO_FREQ, ds_test: xr.Dataset
    ):
        """Get the reference climatology dataset for the variable and season.

        If the reference climatatology does not exist or could not be found, it
        will be considered a model-only run. For this case the test dataset
        is returned as a default value and subsequent metrics calculations will
        only be performed on the original test dataset.

        Parameters
        ----------
        var_key : str
            The key of the variable.
        season : CLIMO_FREQ
            The climatology frequency.
        ds_test : xr.Dataset
            The test dataset, which is returned if the reference climatology
            does not exist or could not be found.

        Returns
        -------
        xr.Dataset
            The reference climatology if it exists or a copy of the test dataset
            if it does not exist.

        Raises
        ------
        RuntimeError
            If `self.data_type` is not "ref".
        """
        # TODO: This logic was carried over from legacy implementation. It
        # can probably be improved on by setting `ds_ref = None` and not
        # performing unnecessary operations on `ds_ref` for model-only runs,
        # since it is the same as `ds_test``.
        if self.data_type == "ref":
            try:
                ds_ref = self.get_climo_dataset(var_key, season)
                self.model_only = False
            except (RuntimeError, IOError):
                ds_ref = ds_test.copy()
                self.model_only = True

                logger.info("Cannot process reference data, analyzing test data only.")
        else:
            raise RuntimeError(
                "`Dataset._get_ref_dataset` only works with "
                f"`self.data_type == 'ref'`, not {self.data_type}."
            )

        return ds_ref

    def get_climo_dataset(self, var: str, season: CLIMO_FREQ) -> xr.Dataset:
        """Get the dataset containing the climatology variable.

        These variables can either be from the test data or reference data.
        If the variable is already a climatology variable, then get it directly
        from the dataset. If the variable is a time series variable, get the
        variable from the dataset and compute the climatology based on the
        selected frequency.

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
        ValueError
            If unable to determine if the variable is a reference or test
            variable and where to find the variable (climatology or time series
            file).
        """
        self.var = var

        if not isinstance(self.var, str) or self.var == "":
            raise ValueError("The `var` argument is not a valid string.")
        if not isinstance(season, str) or season not in CLIMO_FREQS:
            raise ValueError(
                "The `season` argument is not a valid string. Options include: "
                f"{CLIMO_FREQS}"
            )

        if self.is_climo:
            ds = self._get_climo_dataset(season)
        elif self.is_time_series:
            ds = self.get_time_series_dataset(var)
            ds[self.var] = climo(ds, self.var, season)

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
        ds = xr.open_dataset(filepath, use_cftime=True)

        if self.var in ds.variables:
            pass
        elif self.var in self.derived_vars_map:
            ds = self._get_dataset_with_derived_climo_var(ds)
        else:
            raise IOError(
                f"Variable '{self.var}' was not in the file '{filepath}', nor was "
                "it defined in the derived variables dictionary."
            )

        ds = self._squeeze_time_dim(ds)

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
        # filename if the season is in ["ANN", "DJF", "MAM", "JJA", "SON"].
        if season in ["ANN", "DJF", "MAM", "JJA", "SON"]:
            for file in files_in_dir:
                if file.startswith(filename) and season in file:
                    return os.path.join(root_path, file)

        return None

    def _get_dataset_with_derived_climo_var(self, ds: xr.Dataset) -> xr.Dataset:
        """Get the dataset containing the derived variable (`self.var`).

        Parameters
        ----------
        ds: xr.Dataset
            The climatology dataset, whic should contain the source variables
            for deriving the target variable.

        Returns
        -------
        xr.Dataset
            The dataset with the derived variable.
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
        matching_target_var_map = self._get_matching_climo_src_vars(
            ds, target_var, target_var_map
        )
        # Since there's only one set of vars, we get the first and only set
        # of vars from the derived variable dictionary.
        src_var_keys = list(matching_target_var_map.keys())[0]

        # Get the source variable DataArrays and apply the derivation function.
        # Example:
        #   [xr.DataArray(name="PRECC",...), xr.DataArray(name="PRECL",...)]
        src_vars = []
        for var in src_var_keys:
            src_vars.append(ds[var])

        derivation_func = list(matching_target_var_map.values())[0]
        derived_var: xr.DataArray = derivation_func(*src_vars)

        # Add the derived variable to the final xr.Dataset object and return it.
        ds_final = ds.copy()
        ds_final[target_var] = derived_var

        return ds_final

    def _get_matching_climo_src_vars(
        self,
        dataset: xr.Dataset,
        target_var: str,
        target_variable_map: DerivedVariableMap,
    ) -> Dict[Tuple[str, ...], Callable]:
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
        DerivedVariableMap
            The matching dictionary with the key being the source variables
            and the value being the derivation function.

        Raises
        ------
        IOError
            If the datasets for the target variable and source variables were
            not found in the data directory.
        """
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

            if vars_in_file.issuperset(tuple(var_list)):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {tuple(var_list): target_variable_map[var_tuple]}

        raise IOError(
            f"The dataset file has no matching souce variables for {target_var}"
        )

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
            If True, center the time coordinates using time bounds.

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
            ds = xc.center_times(ds)

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
        src_var_keys = list(matching_target_var_map.keys())[0]

        # Unlike the climatology dataset, the source variables for
        # time series data can be found in multiple datasets so a single
        # xr.Dataset object is returned containing all of them.
        ds = self._get_dataset_with_source_vars(src_var_keys)

        # Get the source variable DataArrays.
        # Example:
        #   [xr.DataArray(name="PRECC",...), xr.DataArray(name="PRECL",...)]
        src_vars = [ds[var] for var in src_var_keys]

        # Using the source variables, apply the matching derivation function.
        derivation_func = list(matching_target_var_map.values())[0]
        derived_var: xr.DataArray = derivation_func(*src_vars)

        # Add the derived variable to the final xr.Dataset object and return it.
        ds[target_var] = derived_var

        return ds

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
            if all(self._get_timeseries_filepath(path, var) for var in tuple_of_vars):
                # All of the variables (list_of_vars) have files in data_path.
                # Return the corresponding dict.
                return {tuple_of_vars: target_var_map[tuple_of_vars]}

        # None of the entries in the derived variables dictionary are valid,
        # so try to get the dataset for the variable directly.
        # Example file name: {var}_{start_yr}01_{end_yr}12.nc.
        if self._get_timeseries_filepath(path, self.var):
            return {(self.var,): lambda x: x}

        raise IOError(
            f"Neither does {self.var} nor the variables in {possible_vars} "
            f"have valid files in {path}."
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
        filename = self._get_timeseries_filepath(self.root_path, var)

        if filename == "":
            raise IOError(
                f"No time series `.nc` file was found for '{var}' in '{self.root_path}'"
            )

        time_slice = self._get_time_slice(filename)

        ds = xr.open_dataset(filename, decode_times=True, use_cftime=True)
        ds_subset = ds.sel(time=time_slice).squeeze()

        return ds_subset

    def _get_timeseries_filepath(self, root_path: str, var_key: str) -> str:
        """Get the matching variable time series filepath.

        This method globs the specified path for all `*.nc` files and attempts
        to find a matching time series filepath for the specified variable.

        Example matching filenames.
            - {var}_{start_yr}01_{end_yr}12.nc
            - {self.parameters.ref_name}/{var}_{start_yr}01_{end_yr}12.nc

        If there are multiple files that exist for a variable (with different
        start_yr or end_yr), return an empty string ("").

        Parameters
        ----------
        root_path : str
            The root path containing `.nc` files. The `.nc` files can be nested
            in sub-directories within the root path.
        var_key : str
            The variable key used to find the time series file.

        Returns
        -------
        str
            The variable's time series filepath if a match is found. If
            a match is not found, an empty string ("") is returned.

        Raises
        ------
        IOError
            Multiple time series files found for the specified variable.
        IOError
            Multiple time series files found for the specified variable.
        """
        # The filename pattern for matching using regex.
        if self.parameter.sets[0] in ["arm_diags"]:
            # Example: "ts_global_200001_200112.nc"
            site = getattr(self.parameter, "regions", "")
            filename_pattern = var_key + "_" + site[0] + TS_EXT_FILEPATTERN
        else:
            # Example: "ts_200001_200112.nc"
            filename_pattern = var_key + TS_EXT_FILEPATTERN

        # Attempt 1 -  try to find the file directly in `data_path`
        # Example: {path}/ts_200001_200112.nc"
        match = self._get_matching_time_series_filepath(
            root_path, var_key, filename_pattern
        )

        # Attempt 2 -  try to find the file in the `ref_name` directory, which
        # is nested in `data_path`.
        # Example: {path}/*/{ref_name}/*/ts_200001_200112.nc"
        ref_name = getattr(self.parameter, "ref_name", None)
        if match is None and ref_name is not None:
            match = self._get_matching_time_series_filepath(
                root_path, var_key, filename_pattern, ref_name
            )

        # If there are still no matching files, return an empty string.
        if match is None:
            return ""

        return match

    def _get_matching_time_series_filepath(
        self,
        root_path: str,
        var_key: str,
        filename_pattern: str,
        ref_name: str | None = None,
    ) -> str | None:
        """Get the matching filepath.

        Parameters
        ----------
        root_path : str
            The root path containing `.nc` files. The `.nc` files can be nested
            in sub-directories within the root path.
        var_key : str
            The variable key used to find the time series file.
        filename_pattern : str
            The filename pattern (e.g., "ts_200001_200112.nc").
        ref_name : str | None, optional
            The directory name storing reference files, by default None.

        Returns
        -------
        str | None
            The matching filepath if it exists, or None if it doesn't.

        Raises
        ------
        IOError
            If there are more than one matching filepaths for a variable.
        """
        if ref_name is None:
            # Example: {path}/ts_200001_200112.nc"
            glob_path = os.path.join(root_path, "*.*")
            filepath_pattern = os.path.join(glob_path, filename_pattern)
        else:
            # Example: {path}/{ref_name}/ts_200001_200112.nc"
            glob_path = os.path.join(root_path, ref_name, "*.*")
            filepath_pattern = os.path.join(root_path, ref_name, filename_pattern)

        # Sort the filepaths and loop over them, then check if there are any
        # regex matches using the filepath pattern.
        filepaths = sorted(glob.glob(glob_path))
        matches = [f for f in filepaths if re.search(filepath_pattern, f)]

        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            raise IOError(
                (
                    "There are multiple time series files found for the variable "
                    f"'{var_key}' in '{root_path}' but only one is supported. "
                )
            )

        return None

    def _get_time_slice(self, filename: str) -> slice:
        """Get time slice to subset a dataset.

        Parameters
        ----------
        filename : str
            The filename.

        Returns
        -------
        slice
            A slice object with a start and end time in the format "YYYY-MM-DD".

        Raises
        ------
        ValueError
            If invalid date range specified for test/reference time series data.
        """
        start_year = int(self.start_yr)
        end_year = int(self.end_yr)

        if self.is_sub_monthly:
            start_time = f"{start_year}-01-01"
            end_time = f"{str(int(end_year) + 1)}-01-01"
        else:
            start_time = f"{start_year}-01-15"
            end_time = f"{end_year}-12-15"

        # Get the available start and end years from the file name.
        # Example: {var}_{start_yr}01_{end_yr}12.nc
        var_start_year = int(filename.split("/")[-1].split("_")[-2][:4])
        var_end_year = int(filename.split("/")[-1].split("_")[-1][:4])

        if start_year < var_start_year:
            raise ValueError(
                "Invalid year range specified for test/reference time series data: "
                f"start_year ({start_year}) < var_start_yr ({var_start_year})."
            )
        elif end_year > var_end_year:
            raise ValueError(
                "Invalid year range specified for test/reference time series data: "
                f"end_year ({end_year}) > var_end_yr ({var_end_year})."
            )

        return slice(start_time, end_time)

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
        try:
            ds_land_frac = self.get_climo_dataset(LAND_FRAC_KEY, season)  # type: ignore
            ds_ocean_frac = self.get_climo_dataset(OCEAN_FRAC_KEY, season)  # type: ignore
        except IOError as e:
            logger.info(
                f"{e}. Using default land sea mask located at `{LAND_OCEAN_MASK_PATH}`."
            )

            ds_mask = xr.open_dataset(LAND_OCEAN_MASK_PATH)
            ds_mask = self._squeeze_time_dim(ds_mask)
        else:
            ds_mask = xr.merge([ds_land_frac, ds_ocean_frac])

        return ds_mask

    def _squeeze_time_dim(self, ds: xr.Dataset) -> xr.Dataset:
        """Squeeze single coordinate climatology time dimensions.

        For example, "ANN" averages over the year and collapses the time dim.
        Parameters
        ----------
        ds : xr.Dataset
            _description_

        Returns
        -------
        xr.Dataset
            _description_
        """
        dim = xc.get_dim_keys(ds[self.var], axis="T")
        ds = ds.squeeze(dim=dim)
        ds = ds.drop_vars(dim)

        return ds
