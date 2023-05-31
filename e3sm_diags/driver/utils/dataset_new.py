"""
Get a variable from input data (either reference or test data).
This data can either be climatology files or timeseries files.
Derived variables are also supported.
"""
import collections
import fnmatch
import glob
import os
import re
from typing import Callable, Dict, Literal, Tuple

import cdms2
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.acme_new import (
    DERIVED_VARIABLES,
    DerivedVariableMap,
    DerivedVariablesMap,
)
from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQ, CLIMO_FREQS, climo
from e3sm_diags.parameter.core_parameter import CoreParameter


class Dataset:
    def __init__(
        self,
        parameter: CoreParameter,
        type: Literal["ref", "test"],
    ):
        # The CoreParameter object with a list of parameters.
        self.parameter = parameter

        # The type of data for the Dataset object to store.
        self.type = type

        # The path, start year, and end year based on the dataset type.
        if self.type == "ref":
            self.data_path = self.parameter.reference_data_path
        elif self.type == "test":
            self.data_path = self.parameter.test_data_path
        else:
            raise ValueError(
                f"The `type` ({self.type}) for this Dataset object is invalid."
                "Valid options include 'ref' or 'test'."
            )

        # The start and end year attributes is different for the
        # area_mean_time_series parameter.
        if self.parameter.sets[0] in ["area_mean_time_series"]:
            self.start_yr = self.parameter.start_yr  # type: ignore
            self.end_yr = self.parameter.end_yr  # type: ignore
        elif self.type == "ref":
            self.start_yr = self.parameter.ref_start_yr  # type: ignore
            self.end_yr = self.parameter.ref_end_yr  # type: ignore
        elif self.type == "test":
            self.start_yr = self.parameter.test_start_yr  # type: ignore
            self.end_yr = self.parameter.test_end_yr  # type: ignore

        # The derived variables defined in E3SM Diags. If the CoreParameter
        # object contains additional user derived variables, they are added
        # to `self.derived_vars`.
        self.derived_vars_map = self._get_derived_vars_map()

        # Whether the data is sub-monthly or not.
        self.is_sub_monthly = False
        if self.parameter.sets[0] in ["diurnal_cycle", "arm_diags"]:
            self.is_sub_monthly = True

    @property
    def is_time_series(self):
        if self.type == "ref":
            return self.parameter.ref_timeseries_input
        elif self.type == "test":
            return self.parameter.test_timeseries_input

    @property
    def is_climo(self):
        return not self.is_time_series

    def get_name_and_yrs(self, season: str = ""):
        name = self._get_name()
        yrs_averaged = self._get_yrs(season)
        if yrs_averaged:
            name_yrs = "{} ({})".format(name, yrs_averaged)
        else:
            name_yrs = name

        return name_yrs

    def _get_name(self):
        if self.type == "test":
            if self.parameter.short_test_name:
                name = self.parameter.short_test_name
            else:
                name = self.parameter.test_name
        elif self.type == "ref":
            if self.parameter.short_ref_name:
                name = self.parameter.short_ref_name
            elif self.parameter.reference_name != "":
                # parameter.ref_name is used to search though the reference data directories.
                # parameter.reference_name is printed above ref plots.
                name = self.parameter.reference_name
            else:
                name = self.parameter.ref_name
        return name

    def _get_yrs(self, season=""):
        if self.is_climo:
            try:
                yrs_averaged = self._get_attr_from_climo("yrs_averaged", season)
            except Exception:
                yrs_averaged = ""
        elif self.is_time_series:
            yrs_averaged = "{}-{}".format(self.start_yr, self.end_yr)

        return yrs_averaged

    def _get_derived_vars_map(self) -> DerivedVariablesMap:
        """Get the defined derived variables.

        If the user-defined derived variables is in the input parameters, append
        parameters.derived_variables to the correct part of the derived
        variables dictionary.

        Returns
        -------
        DerivedVariablesMap
            A dictionary mapping the key of a derived variable to an ordered
            dictionary that maps a tuple of source variable(s) to a derivation
            function.
        """
        dvars = DERIVED_VARIABLES.copy()
        user_dvars = getattr(self.parameter, "derived_variables")

        if user_dvars is not None:
            # If the user-defined derived vars already exist, create a
            # new OrderedDict that combines the user-defined entries with the
            # existing ones in `e3sm_diags`. The user-defined entry should
            # be the highest priority and must be first in the OrderedDict.
            for key, ordered_dict in user_dvars.items():
                if key in dvars.keys():
                    dvars[key] = collections.OrderedDict(**ordered_dict, **dvars[key])
                else:
                    dvars[key] = ordered_dict

        return dvars

    # --------------------------------------------------------------------------
    # Climatology related methods
    # --------------------------------------------------------------------------
    def get_climo_dataset(self, var: str, season: CLIMO_FREQ) -> xr.Dataset:
        """Get climatology variables from climatology datasets.

        These variables can either be from the test data or reference data.
        For a given season, get the variable and any extra variables and run
        the climatology on them.

        If the variable is a climatology variable then get it directly
        from the dataset. If the variable is a time series variable, get the
        variable from the dataset and compute the climatology.

        Parameters
        ----------
        var : str
            The variable name
        season : CLIMO_FREQ, optional
            The season for calculation climatology.

        Returns
        -------
        xr.Dataset
            The dataset containing the climatology of the variable.

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
            filename = self._get_climo_filename(season)
            ds = self._get_climo_dataset(filename)
        elif self.is_time_series:
            ds = self.get_time_series_dataset(var)
            ds[self.var] = climo(ds, self.var, season)

        return ds

    def _get_climo_filename(self, season: str) -> str:
        """Return the path to the climatology file

        Parameters
        ----------
        season : str
            The season.

        Returns
        -------
        str
            The path to the climatology file.

        Raises
        ------
        IOError
            If the reference data filepath could not be found.
        IOError
            If the test data filepath could not be found.
        """
        filename = None

        if self.type == "ref":
            data_name = self.parameter.ref_name

            if self.parameter.ref_file != "":
                filename = os.path.join(self.data_path, self.parameter.ref_file)

        elif self.type == "test":
            data_name = self.parameter.test_name

            if hasattr(self.parameter, "test_file"):
                filename = os.path.join(self.data_path, self.parameter.test_file)

        if filename is None:
            filename = self._find_climo_filename(self.data_path, data_name, season)
        elif not os.path.exists(filename):
            raise IOError("File not found: {}".format(filename))

        return filename

    def _find_climo_filename(self, path: str, data_name: str, season: str):
        """
        For climo files, return the path of the file based on the parameters.
        If the file isn't found, try looking for it in path/data_name/ dir as well.
        """
        # First attempt: try to find the climatology file based on season.
        # (e.g., {path}/{data_name}_{season})
        filename = self._find_climo_filename_with_season(path, data_name, season)

        # Second attempt: try looking for the file nested in a folder, based on
        # the test_name (e.g., {path}/{data_name})
        if filename is None:
            pth = os.path.join(path, data_name)

            if os.path.exists(pth):
                filename = self._find_climo_filename_with_season(pth, data_name, season)

        # If absolutely no filename was found, then raise an error.
        if filename is None:
            raise IOError(f"No file found for {data_name} and {season} in {path}")

        return filename

    def _find_climo_filename_with_season(
        self, path_name: str, data_name: str, season: str
    ):
        """
        Locate climatology file name based on data_name and season.
        """
        dir_files = sorted(os.listdir(path_name))
        for filename in dir_files:
            if filename.startswith(data_name + "_" + season):
                return os.path.join(path_name, filename)

        # The below is only ran on model data, because a shorter name is passed
        # into this software. Won't work when use month name such as '01' as
        # season.
        for filename in dir_files:
            if season in ["ANN", "DJF", "MAM", "JJA", "SON"]:
                if filename.startswith(data_name) and season in filename:
                    return os.path.join(path_name, filename)

        return None

    def _get_climo_dataset(self, filepath: str):
        """For a given season and climo input data, get the variable (self.var)."""
        ds = xr.open_dataset(filepath)

        if self.var in self.derived_vars_map:
            ds = self._get_dataset_with_derived_climo_var(filepath)
        elif self.var in ds.variables:
            ds[self.var] = ds[self.var].squeeze(axis=1)
        else:
            raise IOError(
                f"Variable '{self.var}' was not in the file {ds.uri}, nor was "
                "it defined in the derived variables dictionary."
            )

        return ds

    def _get_dataset_with_derived_climo_var(self, filepath: str) -> xr.Dataset:
        """Get the dataset containing the derived variable (`self.var`).

        Parameters
        ----------
        filepath: str
            The filepath to the dataset.

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

        # The climatology dataset is opened up directly and should contain
        # the source variables for deriving the target variable.
        ds = xr.open_dataset(filepath)

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

        # Get the source variable DataArrays and squeeze down the singleton
        # time axis since it is a climatology dataset that has been averaged
        # over the time axis.
        # Example:
        #   [xr.DataArray(name="PRECC",...), xr.DataArray(name="PRECL",...)]
        src_vars = [ds[var].squeeze(axis=1) for var in src_var_keys]

        # Using the source variables, apply the matching derivation function.
        derivation_func = list(matching_target_var_map.values())[0]
        derived_var: xr.DataArray = derivation_func(*src_vars)

        # Add the derived variable to the final xr.Dataset object and return it.
        ds[target_var] = derived_var

        return ds

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
        # ex: [('pr',), ('PRECC', 'PRECL')]
        possible_vars = list(target_variable_map.keys())

        # Add support for wild card `?` in variable strings: ex ('bc_a?DDF', 'bc_c?DDF')
        for list_of_vars in possible_vars:
            matched_var_list = list(list_of_vars).copy()
            for var_list in list_of_vars:
                if "?" in var_list:
                    matched_var_list += fnmatch.filter(list(vars_in_file), var_list)
                    matched_var_list.remove(var_list)

            if vars_in_file.issuperset(tuple(matched_var_list)):
                # All of the variables (list_of_vars) are in data_file.
                # Return the corresponding dict.
                return {tuple(matched_var_list): target_variable_map[list_of_vars]}

        # None of the entries in the derived vars dictionary work,
        # so try to get the var directly.
        # Only try this if var actually exists in data_file.
        if target_var in dataset.data_vars.keys():
            # The below will just cause var to get extracted from the data_file.
            return {(target_var,): lambda x: x}

        raise IOError(
            f"Neither does {target_var} nor the variables in {possible_vars} "
            f"exist in the file {dataset.uri}."
        )

    def _get_attr_from_climo(self, attr, season):
        # TODO: Refactor this method.
        """
        For the given season, get the global attribute from the corresponding climo file.
        """
        if self.is_time_series:
            raise TypeError("Cannot get a global attribute from timeseries files.")

        filename = self._get_climo_filename(season)

        with cdms2.open(filename) as f:
            return f.getglobal(attr)

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
            The time series variable.
        single_point : bool, optional
            Single point indicating sub monthly, by default False

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

        # TODO: Consider making ds a class attribute.
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
            self.data_path, target_var_map
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
        filename = self._get_timeseries_filepath(self.data_path, var)

        if filename == "":
            raise IOError(
                f"No time series `.nc` file was found for '{var}' in '{self.data_path}'"
            )

        time_slice = self._get_time_slice(filename)

        ds = xr.open_dataset(filename, decode_times=True, use_cftime=True)
        ds_subset = ds.sel(time=time_slice).squeeze()

        return ds_subset

    def _get_timeseries_filepath(self, path: str, var_key: str) -> str:
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
        path : str
            The path containing `.nc` files.
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
        # NOTE: Everything between '{var}_' and '.nc' in a time series file is
        # always 13 characters.
        if self.parameter.sets[0] in ["arm_diags"]:
            # Example: "ts_global_200001_200112.nc"
            site = getattr(self.parameter, "regions", "")
            filename_pattern = var_key + "_" + site[0] + r"_.{13}.nc"
        else:
            # Example: "ts_200001_200112.nc"
            filename_pattern = var_key + r"_.{13}.nc"

        # The root directory for the data.
        path = os.path.join(path, "*.*")

        # Attempt 1 -  try to find the file directly in `data_path`
        # Example: {path}/ts_200001_200112.nc"
        filepath_pattern = os.path.join(path, filename_pattern)
        filepaths = sorted(glob.glob(path))
        matches = [f for f in filepaths if re.search(filepath_pattern, f)]

        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            raise IOError(
                (
                    "There are multiple time series files found for the variable "
                    f"'{var_key}' in '{path}' but only one is supported. "
                )
            )

        # Attempt 2 -  try to find the `ref_name` directory nested in `data_path`.
        # Example: {path}/{ref_name}/ts_200001_200112.nc"
        ref_name = getattr(self.parameter, "ref_name", "")
        ref_name_path = os.path.join(path, ref_name, "*.*")
        filepaths = sorted(glob.glob(ref_name_path))

        # Get all of the files that match the pattern
        filepath_pattern = os.path.join(ref_name_path, ref_name, filename_pattern)
        matches = [f for f in filepaths if re.search(filepath_pattern, f)]

        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            raise IOError(
                (
                    "There are multiple time series files found for the variable "
                    f"'{var_key}' in '{path}' but only one is supported. "
                )
            )

        return ""

    def _get_time_slice(self, filename: str) -> slice:
        """Get time slice to subset a dataset.


        Parameters
        ----------
        filename : str
            The filename

        Returns
        -------
        slice
            A slice object with a start and end time in the format "YYYY-MM-DD".

        Raises
        ------
        RuntimeError
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
                f"start_year={start_year}>{var_start_year}=var_start_yr"
            )
        elif end_year > var_end_year:
            raise ValueError(
                "Invalid year range specified for test/reference time series data: "
                f"end_year={end_year}>{var_end_year}=var_end_yr"
            )

        return slice(start_time, end_time)