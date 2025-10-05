from __future__ import annotations

import os
import traceback
from typing import TYPE_CHECKING, Literal, get_args

import uxarray as ux
import xarray as xr

from e3sm_diags.derivations.derivations import DERIVED_VARIABLES, FUNC_NEEDS_TARGET_VAR
from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.logger import _setup_child_logger

if TYPE_CHECKING:
    from collections.abc import Callable

    from e3sm_diags.driver.utils.type_annotations import TimeSelection
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

logger = _setup_child_logger(__name__)


class NativeDataset:
    """
    A class for handling native grid datasets using xarray for raw data
    and uxarray for grid-aware operations.

    NOTE: NativeDataset uses composition instead of inheritance to wrap Dataset
    with additional native-grid specific functionalities. It does not inherit
    from Dataset to avoid confusion with existing Dataset methods and prevent
    tight coupling. If needed, we can refactor to inherit from Dataset or
    create a parent abstract class in the future.
    """

    def __init__(
        self,
        parameter: LatLonNativeParameter,
        data_type: Literal["test", "ref"],
    ):
        # The dataset object (test or reference).
        self.dataset = Dataset(parameter, data_type)

        # The uxarray dataset for grid operations.
        self.grid_dataset: ux.Dataset | None = None

    @property
    def dataset_name(self) -> str:
        return "reference" if self.dataset.data_type == "ref" else "test"

    # --------------------------------------------------------------------------
    # Native-dataset related methods
    # --------------------------------------------------------------------------
    def get_native_dataset(
        self,
        var_key: str,
        season: TimeSelection,
        is_time_slice: bool = False,
        allow_missing: bool = False,
    ) -> xr.Dataset | None:
        """Get the climatology dataset for the variable and season for native grid processing.

        This function handles both test and reference datasets. For reference datasets,
        if the data cannot be found and allow_missing=True, it will return None to
        enable model-only runs.

        This function also stores the data file path in the parameter object
        for native grid visualization.

        Parameters
        ----------

        var_key : str
            The key of the variable.
        season : TimeSelection
            The climatology frequency or time slice string.
        is_time_slice : bool, optional
            If True, treat season as a time slice string rather than climatology
            frequency. Default is False.
        allow_missing : bool, optional
            If True, return None when dataset cannot be loaded instead of raising
            an exception. This enables model-only runs when reference data is
            missing. Default is False.

        Returns
        -------
        xr.Dataset | None
            The climatology dataset if it exists, or None if allow_missing=True
            and the dataset cannot be loaded.

        Raises
        ------
        RuntimeError, IOError
            If the dataset cannot be loaded and allow_missing=False.
        """
        try:
            if is_time_slice:
                ds = self._get_full_native_dataset()
                ds = self._apply_time_slice(ds, season)
            else:
                if season in get_args(ClimoFreq):
                    ds = self.dataset.get_climo_dataset(var_key, season)  # type: ignore
                else:
                    raise ValueError(f"Invalid season for climatology: {season}")

            # Store file path in parameter for native grid processing.
            # Note: For climatology case, get_climo_dataset() already handles file
            # path storage.
            if is_time_slice:
                # For time slices, we know the exact file path we used.
                filepath = self.dataset._get_climo_filepath_with_params()

                if filepath:
                    if self.dataset.data_type == "test":
                        self.dataset.parameter.test_data_file_path = filepath
                    elif self.dataset.data_type == "ref":
                        self.dataset.parameter.ref_data_file_path = filepath

            return ds
        except (RuntimeError, IOError) as e:
            if allow_missing:
                logger.info(
                    f"Cannot process {self.dataset.data_type} data: {e}. Using model-only mode."
                )
                return None
            else:
                raise

    def _get_full_native_dataset(self) -> xr.Dataset:
        """Get the full native dataset without any time averaging.

        This function uses the dataset's file path parameters to directly open
        the raw data file for time slicing operations.

        Parameters
        ----------
        dataset : Dataset
            The dataset object (test or reference).
        var_key : str
            The key of the variable.

        Returns
        -------
        xr.Dataset
            The full dataset with all time steps.

        Raises
        ------
        RuntimeError
            If unable to get the full dataset.
        """
        filepath = self.dataset._get_climo_filepath_with_params()

        if filepath is None:
            raise RuntimeError(
                f"Unable to get file path for {self.dataset.data_type} dataset. "
                f"For time slicing, please ensure that "
                f"{'ref_file' if self.dataset.data_type == 'ref' else 'test_file'} parameter is set."
            )

        if not os.path.exists(filepath):
            raise RuntimeError(f"File not found: {filepath}")

        logger.info(f"Opening full native dataset from: {filepath}")

        try:
            # Open the dataset directly without any averaging
            ds = xr.open_dataset(filepath, decode_times=True)
            logger.info(
                f"Successfully opened dataset with time dimension size: {ds.sizes.get('time', 'N/A')}"
            )

            return ds
        except Exception as e:
            raise RuntimeError(f"Failed to open dataset {filepath}: {e}") from e

    def _apply_time_slice(self, ds: xr.Dataset, time_slice: str) -> xr.Dataset:
        """Apply time slice selection to a dataset.

        Parameters
        ----------
        ds : xr.Dataset
            The input dataset with time dimension.
        time_slice : str
            The time slice specification (e.g., "0:10:2", "5:15", "7").

        Returns
        -------
        xr.Dataset
            The dataset with time slice applied.
        """
        # Parse the time slice string
        time_dim = None

        for dim in ds.dims:
            if str(dim).lower() in ["time", "t"]:
                time_dim = dim
                break

        if time_dim is None:
            logger.warning(
                "No time dimension found in dataset. Returning original dataset."
            )
            return ds

        # Parse slice notation
        if ":" in time_slice:
            # Handle slice notation like "0:10:2" or "5:15" or ":10" or "5:" or "::2"
            parts = time_slice.split(":")

            start = int(parts[0]) if parts[0] else None
            end = int(parts[1]) if len(parts) > 1 and parts[1] else None
            step = int(parts[2]) if len(parts) > 2 and parts[2] else None

            # Apply the slice
            ds_sliced = ds.isel({time_dim: slice(start, end, step)})
        else:
            # Single index
            index = int(time_slice)
            ds_sliced = ds.isel({time_dim: index})

        logger.info(
            f"Applied time slice '{time_slice}' to dataset. "
            f"Original time length: {ds.sizes[time_dim]}, "
            f"Sliced time length: {ds_sliced.sizes.get(time_dim, 1)}"
        )

        return ds_sliced

    # --------------------------------------------------------------------------
    # Grid dataset related methods
    # --------------------------------------------------------------------------
    def get_grid_dataset(self) -> ux.Dataset:
        """Open the dataset using uxarray.

        Returns
        -------
        ux.Dataset
            The opened dataset.
        """
        uxds = None

        try:
            if self.dataset.data_type == "test":
                uxds = ux.open_dataset(
                    self.dataset.parameter.test_grid_file,  # type: ignore
                    self.dataset.parameter.test_data_file_path,
                )
            elif self.dataset.data_type == "ref":
                has_ref_grid = (
                    hasattr(self.dataset.parameter, "ref_grid_file")
                    and self.dataset.parameter.ref_grid_file is not None
                )

                if not has_ref_grid:
                    logger.info(
                        "No ref_grid_file specified. Skipping reference grid loading."
                    )
                else:
                    grid_file = self.dataset.parameter.ref_grid_file  # type: ignore

                # Use ref_data_file_path if available, otherwise use ds_ref
                if (
                    hasattr(self.dataset.parameter, "ref_data_file_path")
                    and self.dataset.parameter.ref_data_file_path
                ):
                    data_source = self.dataset.parameter.ref_data_file_path
                else:
                    data_source = uxds  # type: ignore

                uxds = ux.open_dataset(grid_file, data_source)
        except Exception as e:
            logger.error(f"Failed to load {self.dataset.data_type} native grid: {e}")

            logger.debug(traceback.format_exc())

        self.grid_dataset = uxds

        return uxds

    def _process_variable_derivations(self, var_key: str) -> bool:
        """Process variable derivations following dataset_xr approach."""
        name_suffix = f" in {self.dataset_name} dataset" if self.dataset_name else ""

        # Follow dataset_xr._get_climo_dataset logic:
        # 1. If var is in derived_vars_map, try to derive it
        if var_key in DERIVED_VARIABLES:
            target_var_map = DERIVED_VARIABLES[var_key]
            matching_target_var_map = self._get_matching_src_vars(
                self.grid_dataset, target_var_map
            )

            if matching_target_var_map is not None:
                # Get derivation function and source variables
                derivation_func = list(matching_target_var_map.values())[0]
                src_var_keys = list(matching_target_var_map.keys())[0]

                logger.info(
                    f"Deriving {var_key}{name_suffix} using source variables: {src_var_keys}"
                )

                try:
                    self._apply_derivation_func(
                        self.grid_dataset, derivation_func, src_var_keys, var_key
                    )
                    return True
                except Exception as e:
                    logger.warning(f"Failed to derive {var_key}{name_suffix}: {e}")

        # 2. Check if variable exists directly in dataset
        if var_key in self.grid_dataset.data_vars:  # type: ignore
            return True

        # 3. Variable not found and couldn't be derived
        logger.warning(
            f"Variable {var_key} not found{name_suffix} and could not be derived"
        )
        return False

    def _get_matching_src_vars(self, dataset, target_var_map):
        """Get matching source variables following dataset_xr pattern."""
        for src_vars, func in target_var_map.items():
            if all(v in dataset for v in src_vars):
                return {src_vars: func}

        return None

    def _apply_derivation_func(
        self,
        dataset: ux.Dataset,
        func: Callable,
        src_var_keys: list[str],
        target_var_key: str,
    ):
        """Apply derivation function following dataset_xr pattern."""
        func_args = [dataset[var] for var in src_var_keys]

        if func in FUNC_NEEDS_TARGET_VAR:
            func_args = [target_var_key] + func_args

        derived_var = func(*func_args)
        dataset[target_var_key] = derived_var

        return dataset
