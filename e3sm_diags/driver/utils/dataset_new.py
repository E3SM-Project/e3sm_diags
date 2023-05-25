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
from typing import Any, Callable, Dict, List, Literal, Tuple

import cdms2
import xarray as xr
import xcdat as xc

from e3sm_diags.derivations.acme_new import derived_variables
from e3sm_diags.driver.utils.climo_xr import CLIMO_FREQ, climo
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

        # The path to to the data based on the type of data.
        if self.type == "ref":
            self.data_path = self.parameter.reference_data_path
            self.start_yr = getattr(self.parameter, "ref_start_yr")
            self.end_yr = getattr(self.parameter, "ref_end_yr")
        elif self.type == "test":
            self.data_path = self.parameter.test_data_path
            self.start_yr = getattr(self.parameter, "test_start_yr")
            self.end_yr = getattr(self.parameter, "test_end_yr")
        else:
            raise ValueError(
                f"The `type` ({self.type})for this Dataset object is invalid."
                "Valid options include 'ref' or 'test'."
            )

        # The derived variables defined in E3SM Diags. If the CoreParameter
        # object contains additional user derived variables, they are added
        # to `self.derived_vars`.
        self.derived_vars = self._get_derived_vars()

        # The start and end year attributes is different for the
        # area_mean_time_series parameter.
        if self.parameter.sets[0] in ["area_mean_time_series"]:
            self.start_yr = getattr(self.parameter, "start_yr")
            self.end_yr = getattr(self.parameter, "end_yr")

        # Whether the data is sub-monthly or not.
        self.sub_monthly = False
        if self.parameter.sets[0] in ["diurnal_cycle", "arm_diags"]:
            self.sub_monthly = True

    @property
    def is_time_series(self):
        if self.type == "ref":
            return getattr(self.parameter, "ref_timeseries_input", False)
        elif self.type == "test":
            return getattr(self.parameter, "test_timeseries_input", False)

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

    def _get_derived_vars(self) -> Dict[str, Dict[Any, Any]]:
        """Get the defined derived variables.

        If the user-defined derived variables is in the input parameters, append
        parameters.derived_variables to the correct part of the derived
        variables dictionary.

        Returns
        -------
        # TODO: Update return type
        Dict[str, OrderedDict]
            A dictionary, with the key being the variable name and the value
            being an OrderedDict mapping variables to their function for
            deriving them using other variables.
        """
        derived_vars = derived_variables

        if hasattr(self.parameter, "derived_variables"):
            key_val_pairs = self.parameter.derived_variables.items()
            for derived_var, original_vars in list(key_val_pairs):
                # Append the user-defined vars to the already defined ones.
                if derived_var in self.derived_vars:
                    # Put user-defined derived vars first in the OrderedDict.
                    # That's why we create a new one.
                    new_dict = collections.OrderedDict(original_vars)  # type: ignore
                    # Add all of the default derived vars to the end of new_dict.
                    for k in self.derived_vars[derived_var]:
                        # Don't overwrite the user-defined var with a default derived var.
                        if k in new_dict:
                            continue
                        new_dict[k] = derived_vars[derived_var][k]  # type: ignore
                    derived_vars[derived_var] = new_dict
                # Otherwise, this is a new derived var, so add it as a new entry.
                else:
                    derived_vars[derived_var] = original_vars

        return derived_vars  # type: ignore

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
        RuntimeError
            If the specified variable is invalid or not set.
        RuntimeError
            If the specified season is invalid or not set.
        RuntimeError
            If unable to determine if the variable is a reference or test
            variable and where to find the variable (climatology or time series
            file).
        """
        self.var = var

        if not self.var:
            raise RuntimeError("Variable is invalid.")
        if not season:
            raise RuntimeError("Season is invalid.")

        if self.is_climo:
            filename = self.get_climo_filename(season)
            ds = self._get_climo_dataset(filename)
        elif self.is_time_series:
            ds = self._get_time_series_dataset(self.data_path)
            ds[self.var] = climo(ds, self.var, season)
        else:
            raise RuntimeError(
                "Error when determining what kind (ref or test) variable to "
                "get and where to get it from (climo or timeseries files)."
            )

        return ds

    def get_climo_filename(self, season: str) -> str:
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
        if self.type == "ref":
            data_name = self.parameter.ref_name

            if (
                hasattr(self.parameter, "ref_file")
                and getattr(self.parameter, "ref_file") != ""
            ):
                fnm = os.path.join(self.data_path, self.parameter.ref_file)

                if not os.path.exists(fnm):
                    raise IOError("File not found: {}".format(fnm))

                return fnm
        elif self.type == "test":
            data_name = self.parameter.test_name

            if hasattr(self.parameter, "test_file"):
                fnm = os.path.join(self.data_path, self.parameter.test_file)

                if not os.path.exists(fnm):
                    raise IOError("File not found: {}".format(fnm))

                return fnm

        return self._get_climo_filename(self.data_path, data_name, season)

    def _get_climo_filename(self, path, data_name, season):
        """
        For climo files, return the path of the file based on the parameters.
        If the file isn't found, try looking for it in path/data_name/ dir as well.
        """
        fnm = self._find_climo_file(path, data_name, season)
        if not os.path.exists(fnm):
            # Try looking for the file nested in a folder, based on the test_name.
            pth = os.path.join(path, data_name)
            if os.path.exists(pth):
                fnm = self._find_climo_file(pth, data_name, season)

        if not os.path.exists(fnm):
            raise IOError(f"No file found for {data_name} and {season} in {path}")

        return fnm

    def _find_climo_file(self, path_name, data_name, season):
        """
        Locate climatology file name based on data_name and season.
        """
        dir_files = sorted(os.listdir(path_name))
        for filename in dir_files:
            if filename.startswith(data_name + "_" + season):
                return os.path.join(path_name, filename)
        # The below is only ran on model data, because a shorter name is passed into this software. Won't work when use month name such as '01' as season.
        for filename in dir_files:
            if season in ["ANN", "DJF", "MAM", "JJA", "SON"]:
                if filename.startswith(data_name) and season in filename:
                    return os.path.join(path_name, filename)
        # No file found.
        return ""

    def _get_climo_dataset(self, path: str):
        """For a given season and climo input data, get the variable (self.var)."""
        ds = xr.open_dataset(path)

        if self.var in self.derived_vars:
            derived_var = self._get_derived_dataset(ds, type="climo")
            # Add the derived variable to the dataset for downstream
            # processing (e.g., climatology).
            ds[derived_var.name] = derived_var
        elif self.var in ds.variables:
            ds[self.var] = ds[self.var].squeeze(axis=1)
        else:
            raise RuntimeError(
                f"Variable '{self.var}' was not in the file {ds.uri}, nor was "
                "it defined in the derived variables dictionary."
            )

        return ds

    def _get_attr_from_climo(self, attr, season):
        # TODO: Refactor this method.
        """
        For the given season, get the global attribute from the corresponding climo file.
        """
        if self.is_time_series:
            raise RuntimeError("Cannot get a global attribute from timeseries files.")

        filename = self.get_climo_filename(season)

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
        RuntimeError
            If the dataset is not a time series.
        """
        self.var = var

        if not self.is_time_series:
            raise RuntimeError("You can only use this function with time series data.")

        ds = self._get_time_series_dataset(self.data_path)

        if single_point:
            ds = xc.center_times(ds)

        return ds

    def _get_time_series_dataset(self, path: str) -> xr.Dataset:
        """Get the time series dataset.

        The time series dataset can either include the derived variable from
        other variables, or the matching time series dataset.

        Parameters
        ----------
        path : str
            The path to the dataset.

        Returns
        -------
        xr.Dataset
            The time series dataset.

        Raises
        ------
        RuntimeError
            If the variable(s) don't have a matching file or is not found in the
            derived variables dictionary.
        """
        # Can't iterate through self.var and self.extra_vars as we do in _get_climo_var()
        # b/c the extra_vars must be taken from the same timeseries file as self.var.
        # So once we got a working vars_to_func_dict, we need to use this to get the extra_vars.

        if self.var in self.derived_vars:
            ds = self._get_derived_dataset(path, type="time_series")
            return ds

        try:
            ds = self._get_time_series_dataset_obj(path, self.var)
            return ds
        except RuntimeError:
            pass

        raise RuntimeError(
            f"Variable '{self.var}' doesn't have a file in the directory "
            f"'{path}', nor was it defined in the derived variables "
            "dictionary."
        )

    def _get_time_series_dataset_obj(self, path: str, var_key: str) -> xr.Dataset:
        """Get the time series dataset for a variable.

        This method also parses the start and end time from the dataset filename
        to subset the dataset.

        Parameters
        ----------
        path : str
            The path to the variable's dataset file.
        var_key : str
            The key of the variable.

        Returns
        -------
        xr.Dataset
            The dataset for the variable.
        """
        filename = self._get_timeseries_filepath(path, var_key)
        time_slice = self._get_time_slice(filename)

        ds = xc.open_dataset(filename, add_bounds=False, decode_times=True)
        ds = ds.sel(time=time_slice).squeeze()

        return ds

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
            The variable's time series filepath if it exists.

        Raises
        ------
        RuntimeError
            Multiple time series files found for the specified variable.
        RuntimeError
            Multiple time series files found for the specified variable.
        """
        # TODO: Refactor this method (repeating lines of code).
        # Get all of the nc file paths in data_path.
        path = os.path.join(path, "*.*")
        files = sorted(glob.glob(path))

        # Both .nc and .xml files are supported
        file_fmt = ""
        if len(files) > 0:
            file_fmt = files[0].split(".")[-1]

        # Everything between '{var}_' and '.nc' in a
        # time-series file is always 13 characters.
        if self.parameter.sets[0] in ["arm_diags"]:
            site = getattr(self.parameter, "regions", "")
            re_str = var_key + "_" + site[0] + r"_.{13}." + file_fmt
        else:
            re_str = var_key + r"_.{13}." + file_fmt
        re_str = os.path.join(path, re_str)
        matches = [f for f in files if re.search(re_str, f)]

        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            msg = "For the variable {} you have two timeseries files in the ".format(
                var_key
            )
            msg += "directory: {} This currently isn't supported.".format(path)
            raise RuntimeError(msg)

        # If nothing was found, try looking for the file with
        # the ref_name prepended to it.
        ref_name = getattr(self.parameter, "ref_name", "")
        # path = os.path.join(data_path, ref_name, '*.nc')
        path = os.path.join(path, ref_name, "*.*")
        files = sorted(glob.glob(path))
        # Both .nc and .xml files are supported
        file_fmt = ""
        if len(files) > 0:
            file_fmt = files[0].split(".")[-1]

        # Everything between '{var}_' and '.nc' in a
        # time-series file is always 13 characters.
        re_str = var_key + r"_.{13}." + file_fmt
        re_str = os.path.join(path, ref_name, re_str)
        matches = [f for f in files if re.search(re_str, f)]
        # Again, there should only be one file per var in this new location.
        if len(matches) == 1:
            return matches[0]
        elif len(matches) >= 2:
            msg = "For the variable {} you have two timeseries files in the ".format(
                var_key
            )
            msg += "directory: {} This currently isn't supported.".format(path)
            raise RuntimeError(msg)
        else:
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
        start_year = self.start_yr
        end_year = self.end_yr

        if self.sub_monthly:
            start_time = "{}-01-01".format(start_year)
            end_time = "{}-01-01".format(str(int(end_year) + 1))
        else:
            start_time = "{}-01-15".format(start_year)
            end_time = "{}-12-15".format(end_year)

        # get available start and end years from file name:
        # {var}_{start_yr}01_{end_yr}12.nc
        start_year = int(start_year)
        end_year = int(end_year)
        var_start_year = int(filename.split("/")[-1].split("_")[-2][:4])
        var_end_year = int(filename.split("/")[-1].split("_")[-1][:4])

        if start_year < var_start_year:
            msg = "Invalid year range specified for test/reference time series data: start_year={}<{}=var_start_yr".format(
                start_year, var_start_year
            )
            raise RuntimeError(msg)
        elif end_year > var_end_year:
            msg = "Invalid year range specified for test/reference time series data: end_year={}>{}=var_end_yr".format(
                end_year, var_end_year
            )
            raise RuntimeError(msg)

        return slice(start_time, end_time)

    # --------------------------------------------------------------------------
    # Derived variable related methods
    # --------------------------------------------------------------------------
    def _get_derived_dataset(self, path: str, type: Literal["climo", "time_series"]):
        # Ex: {('PRECC', 'PRECL'): func, ('pr',): func1, ...}, is an OrderedDict.
        possible_vars_and_funcs = self.derived_vars[self.var]

        # Get the first valid variables and functions from possible vars.
        # Ex: {('PRECC', 'PRECL'): func}
        # These are checked to be in data_file.
        if type == "climo":
            ds = xr.open_dataset(path)
            vars_to_func_dict = self._get_first_valid_vars_climo(
                possible_vars_and_funcs, path, self.var
            )
            # Given a dictionary in the form {(vars): func}, get the vars from the
            # xr.Dataset as xarray.DataArray objects. Since there's only one set of
            # vars, we get the first and only set of vars from the dictionary.
            vars_to_get = list(vars_to_func_dict.keys())[0]
            variables = [ds[var].squeeze(axis=1) for var in vars_to_get]

        elif type == "time_series":
            vars_to_func_dict = self._get_first_valid_vars_timeseries(
                possible_vars_and_funcs, self.var
            )
            # Open the files of the variables and get the xarray.DataArrays.
            # Ex: [PRECC, PRECL], where both are DataArrays.
            ds, variables = self._get_original_vars_timeseries(vars_to_func_dict, path)

        # Derive the variable using the list of variables.
        # Call the first matching function for deriving the variable.
        # Example: {('PRECC', 'PRECL'): func}.
        func = list(vars_to_func_dict.values())[0]
        derived_var: xr.DataArray = func(*variables)

        return derived_var

    def _get_first_valid_vars_climo(self, vars_to_func_dict, data_file, var):
        """
        Given an OrderedDict of a list of variables to a function
             ex: {('PRECC', 'PRECL'): func, ('var2',): func2},
        return the first valid {(vars): func} where the vars are in data_file.

        var is the actual variable the user requested.
        If none of the derived variables work, we try to just get this from the data_file.
        """
        vars_in_file = set(data_file.variables)
        possible_vars = list(
            vars_to_func_dict.keys()
        )  # ex: [('pr',), ('PRECC', 'PRECL')]

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
                return {tuple(matched_var_list): vars_to_func_dict[list_of_vars]}

        # None of the entries in the derived vars dictionary work,
        # so try to get the var directly.
        # Only try this if var actually exists in data_file.
        if var in data_file.variables:
            # The below will just cause var to get extracted from the data_file.
            return {(var,): lambda x: x}

        raise RuntimeError(
            f"Neither does {var} nor the variables in {possible_vars} "
            f"exist in the file {data_file.uri}."
        )

    def _get_first_valid_vars_timeseries(
        self, vars_to_func_dict: Dict[List[str], Callable], path: str
    ) -> Dict[Any, Any]:
        """
        Given an OrderedDict of a list of variables to a function
            ex: {('PRECC', 'PRECL'): func, ('var2',): func2},
        return the first valid {(vars): func} where the vars are variables from files in the form:
            {var}_{start_yr}01_{end_yr}12.nc
        located in data_path.

        If none of the derived variables work, we try to just get self.var in a file like:
            {self.var}_{start_yr}01_{end_yr}12.nc
        located in data_path.
        """
        possible_vars = list(
            vars_to_func_dict.keys()
        )  # ex: [('pr',), ('PRECC', 'PRECL')]

        for list_of_vars in possible_vars:
            # Check that there are files in data_path that exist for all variables in list_of_vars.
            if all(self._get_timeseries_filepath(path, var) for var in list_of_vars):
                # All of the variables (list_of_vars) have files in data_path.
                # Return the corresponding dict.
                return {list_of_vars: vars_to_func_dict[list_of_vars]}

        # None of the entries in the derived vars dictionary are valid,
        # so try to get the var directly.
        # Only try this if there is a corresponding file for var in data_path.
        if self._get_timeseries_filepath(path, self.var):
            # The below will just cause var to get extracted in {var}_{start_yr}01_{end_yr}12.nc.
            return {(self.var,): lambda x: x}

        # Otherwise, there's no way to get the variable.
        raise RuntimeError(
            f"Neither does {self.var} nor the variables in {possible_vars} "
            f"have valid files in {path}."
        )

    def _get_original_vars_timeseries(
        self, vars_to_func_dict: Dict[List[str], Callable], data_path: str
    ) -> Tuple[xr.Dataset, List[xr.DataArray]]:
        """Get the variables from datasets in the specified path.

        Parameters
        ----------
        vars_to_func_dict : Dict[List[str], Callable]
            A dictionary mapping derived variables to their functions. These
            variables are validated to be in the path through
            `_get_first_valid_vars_timeseries().`
        data_path : str
            The path to the datasets.

        Returns
        -------
        Tuple[xr.Dataest, List[xr.DataArray]]
            The dataset and list of variables.
        """
        # Since there's only one set of vars, we get the first and only set of
        # vars from the derived variables dictionary.
        vars_to_get = list(vars_to_func_dict.keys())[0]

        datasets = []
        variables = []

        for var in vars_to_get:
            ds = self._get_time_series_dataset_obj(data_path, var)
            datasets.append(ds)

            v = ds[var].copy()
            variables.append(v)

        ds = xr.merge(datasets)

        return ds, variables
