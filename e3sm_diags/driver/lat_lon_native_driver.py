from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import uxarray as ux
import xarray as xr

from e3sm_diags.driver.utils.climo_xr import ClimoFreq
from e3sm_diags.driver.utils.dataset_xr import Dataset
from e3sm_diags.driver.utils.regrid import (
    _apply_land_sea_mask,
    _subset_on_region,
    get_z_axis,
    has_z_axis,
    regrid_z_axis_to_plevs,
)
from e3sm_diags.driver.utils.type_annotations import MetricsDict
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.metrics.metrics import spatial_avg, std

logger = _setup_child_logger(__name__)

if TYPE_CHECKING:
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

# The default value for metrics if it is not calculated. This value was
# preserved from the legacy CDAT codebase because the viewer expects this
# value for metrics that aren't calculated.
# TODO: Update `lat_lon_viewer.py` to handle missing metrics with None value.
METRICS_DEFAULT_VALUE = 999.999


def run_diag(parameter: LatLonNativeParameter) -> LatLonNativeParameter:
    """Get metrics for the lat_lon_native diagnostic set.

    This function loops over each variable, season, pressure level (if 3-D),
    and region.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter for the diagnostic.

    Returns
    -------
    LatLonNativeParameter
        The parameter for the diagnostic with the result (completed or failed).

    Raises
    ------
    RuntimeError
        If the dimensions of the test and reference datasets are not aligned
        (e.g., one is 2-D and the other is 3-D).
    """
    variables = parameter.variables
    seasons = parameter.seasons
    ref_name = getattr(parameter, "ref_name", "")
    regions = parameter.regions

    # Variables storing xarray `Dataset` objects start with `ds_` and
    # variables storing e3sm_diags `Dataset` objects end with `_ds`. This
    # is to help distinguish both objects from each other.
    test_ds = Dataset(parameter, data_type="test")
    ref_ds = Dataset(parameter, data_type="ref")

    for var_key in variables:
        logger.info("Variable: {}".format(var_key))
        parameter.var_id = var_key

        for season in seasons:
            parameter._set_name_yrs_attrs(test_ds, ref_ds, season)

            ds_test = test_ds.get_climo_dataset(var_key, season)
            ds_ref = _get_ref_climo_dataset(ref_ds, var_key, season)
            ds_land_sea_mask: xr.Dataset = test_ds._get_land_sea_mask(season)

            # Debug information about the test dataset
            logger.info("Test dataset info:")
            logger.info(f"Variables: {list(ds_test.variables)}")
            if hasattr(ds_test, 'file_path'):
                logger.info(f"Dataset file path: {ds_test.file_path}")
            if hasattr(ds_test, 'filepath'):
                logger.info(f"Dataset filepath: {ds_test.filepath}")
            if hasattr(ds_test, '_file_obj') and hasattr(ds_test._file_obj, 'name'):
                logger.info(f"Dataset file object name: {ds_test._file_obj.name}")
            try:
                for name, var in ds_test.variables.items():
                    if hasattr(var, 'file'):
                        logger.info(f"Variable {name} file: {var.file}")
                        break
            except:
                pass

            # Load the native grid information for test data
            uxds_test = None
            if parameter.test_grid_file:
                try:
                    logger.info(f"Loading test native grid from: {parameter.test_grid_file}")
                    # When loading the dataset, include the test data to map it onto the grid
                    uxds_test = ux.open_dataset(parameter.test_grid_file, parameter.test_data_file_path)
                    logger.info("Successfully loaded test native grid data with uxarray")

                    # Verify variables in uxarray dataset
                    if var_key not in uxds_test:
                        logger.warning(f"Variable {var_key} not found in test uxarray dataset!")
                        logger.info(f"Available variables: {list(uxds_test.data_vars)}")

                        # Try to find the variable with a different name or naming convention
                        # For example, some files might use PRECT while others use pr
                        possible_matches = []
                        for data_var in uxds_test.data_vars:
                            if var_key.lower() in data_var.lower() or data_var.lower() in var_key.lower():
                                possible_matches.append(data_var)

                        if possible_matches:
                            logger.info(f"Possible variable matches: {possible_matches}")

                except Exception as e:
                    logger.error(f"Failed to load test native grid: {e}")
                    import traceback
                    logger.error(traceback.format_exc())
                    uxds_test = None

            # Load the native grid information for reference data
            uxds_ref = None
            if ds_ref is not None and parameter.ref_grid_file:
                try:
                    logger.info(f"Loading reference native grid from: {parameter.ref_grid_file}")
                    # When loading the dataset, include the reference data to map it onto the grid
                    uxds_ref = ux.open_dataset(parameter.ref_grid_file, ds_ref)
                    logger.info("Successfully loaded reference native grid data with uxarray")

                    # Verify variables in uxarray dataset
                    if var_key not in uxds_ref:
                        logger.warning(f"Variable {var_key} not found in reference uxarray dataset!")
                        logger.info(f"Available variables: {list(uxds_ref.data_vars)}")
                except Exception as e:
                    logger.error(f"Failed to load reference native grid: {e}")
                    import traceback
                    logger.error(traceback.format_exc())
                    uxds_ref = None

            if ds_ref is None:
                is_vars_3d = has_z_axis(ds_test[var_key])

                if not is_vars_3d:
                    _run_diags_2d_model_only(
                        parameter,
                        ds_test,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                    )
                else:
                    _run_diags_3d_model_only(
                        parameter,
                        ds_test,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                    )
            else:
                is_vars_3d, is_dims_diff = _check_var_dims(ds_test, ds_ref, var_key)

                if is_dims_diff:
                    raise RuntimeError(
                        "Dimensions of the two variables are different. Aborting."
                    )
                elif not is_vars_3d:
                    _run_diags_2d(
                        parameter,
                        ds_test,
                        ds_ref,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                        uxds_ref,  # Use the consistent naming
                    )
                elif is_vars_3d:
                    _run_diags_3d(
                        parameter,
                        ds_test,
                        ds_ref,
                        ds_land_sea_mask,
                        season,
                        regions,
                        var_key,
                        ref_name,
                        uxds_test,  # Use the consistent naming
                        uxds_ref,  # Use the consistent naming
                    )

    return parameter


def _save_plot(
    parameter: LatLonNativeParameter,
    var_key: str,
    region: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Save plots using native grid datasets directly.

    This simplified function creates plots using uxarray dataset(s) without needing
    additional dataset parameters.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    var_key : str
        The variable key.
    region : str
        The region name, used to determine map extents.
    uxds_test : ux.dataset.UxDataset, optional
        The test native grid dataset.
    uxds_ref : ux.dataset.UxDataset, optional
        The reference native grid dataset.
    """
    import os

    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    import matplotlib.colors as mcolors
    import matplotlib.pyplot as plt

    # Import required modules
    from e3sm_diags.derivations.default_regions_xr import REGION_SPECS

    # Check if uxarray dataset and variable are available for plotting
    if uxds_test is None or var_key not in uxds_test:
        logger.error(f"Cannot plot native grid data. Either uxds_test is None or {var_key} not in dataset")
        if uxds_test is not None:
            logger.error(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Create figure
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Get region information for setting map extents
    region_specs = REGION_SPECS.get(region, None)

    # Determine projection and extents based on region
    projection = ccrs.PlateCarree()

    # Set map bounds based on region
    if region_specs:
        lat_bounds = region_specs.get('lat', (-90, 90))
        lon_bounds = region_specs.get('lon', (0, 360))
        is_global_domain = lat_bounds == (-90, 90) and lon_bounds == (0, 360)
    else:
        lat_bounds = (-90, 90)
        lon_bounds = (0, 360)
        is_global_domain = True

    # Add test data subplot
    logger.info(f"Creating test subplot for {var_key}")
    ax = fig.add_subplot(1, 1, 1, projection=projection)
    ax.set_title(f"{parameter.test_name_yrs}\n{parameter.test_title}")

    try:
        # Debug information
        logger.info(f"Plotting variable {var_key} from uxarray dataset")
        logger.info(f"Variable shape: {uxds_test[var_key].shape}")
        logger.info(f"Variable min/max: {uxds_test[var_key].min().item()}/{uxds_test[var_key].max().item()}")

        # Get data range for normalization
        vmin = uxds_test[var_key].min().item()
        vmax = uxds_test[var_key].max().item()

        # Create a normalization for the colormap
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

        # Create a PolyCollection from the native grid data
        logger.info("Creating PolyCollection")
        periodic_elements = "split" if parameter.split_periodic_elements else None

        # Create the PolyCollection
        pc = uxds_test[var_key].squeeze().to_polycollection(
            periodic_elements=periodic_elements,
        )
        # disables grid lines
        pc.set_antialiased(False)

        pc.set_cmap("plasma")

        fig, ax = plt.subplots(
            1,
            1,
            figsize=(10, 5),
            facecolor="w",
            constrained_layout=True,
            subplot_kw=dict(projection=ccrs.PlateCarree()),
        )

        ax.add_feature(cfeature.COASTLINE)
        ax.add_feature(cfeature.BORDERS)

        ax.add_collection(pc)
        ax.set_global()


        ## Explicitly set the array data for colormapping
        #logger.info("Setting up PolyCollection properties")
        #pc.set_array(uxds_test[var_key].values)
        #pc.set_cmap(parameter.test_colormap)
        #pc.set_norm(norm)
        #
        ## Set styling options
        #pc.set_antialiased(parameter.antialiased)
        #pc.set_edgecolor('none')  # Hide grid edges by default
        #
        ## Add the collection
        #logger.info("Adding PolyCollection to axes")
        #ax.add_collection(pc)
        #
        ## Add map features
        #ax.coastlines(linewidth=0.5)
        #ax.add_feature(cfeature.BORDERS, linewidth=0.3)
        #
        ## Set map extent
        #if is_global_domain:
        #    logger.info("Setting global view")
        #    ax.set_global()
        #else:
        #    logger.info(f"Setting map extent to {lon_bounds}, {lat_bounds}")
        #    ax.set_extent([lon_bounds[0], lon_bounds[1], lat_bounds[0], lat_bounds[1]],
        #                 crs=ccrs.PlateCarree())
        #
        ## Add a colorbar
        #logger.info("Adding colorbar")
        #cbar = plt.colorbar(pc, ax=ax, orientation='horizontal',
        #                   pad=0.05, shrink=0.8)
        #cbar.set_label(uxds_test[var_key].attrs.get('units', ''))

    except Exception as e:
        logger.error(f"Error plotting with uxarray: {e}")
        # Print the full exception traceback for debugging
        import traceback
        logger.error(traceback.format_exc())

    # Save the plot - use the same directory structure as in io._get_output_dir
    # This is critical to match the path structure that the viewer expects
    output_dir = os.path.join(parameter.results_dir, parameter.current_set, parameter.case_id)
    os.makedirs(output_dir, exist_ok=True)

    # Log the directory structure for debugging
    logger.info(f"Saving plot to directory: {output_dir}")
    logger.info(f"Output filename base: {parameter.output_file}")

    for fmt in parameter.output_format:
        filename = os.path.join(output_dir, f"{parameter.output_file}.{fmt}")
        fig.savefig(filename)
        logger.info(f"Plot saved to: {filename}")

    plt.close(fig)

def _run_diags_2d_model_only(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
):
    """Run a model-only diagnostics on a 2D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    """
    # Debug the uxarray dataset
    if uxds_test is not None:
        logger.info("===== TEST UXDS DEBUG INFO =====")
        logger.info(f"Type: {type(uxds_test)}")
        logger.info(f"Variables: {list(uxds_test.data_vars)}")
        logger.info(f"Has face attribute: {hasattr(uxds_test, 'face')}")
        if var_key in uxds_test:
            logger.info(f"Variable {var_key} shape: {uxds_test[var_key].shape}")
            logger.info(f"Variable {var_key} min/max: {uxds_test[var_key].min().item()}/{uxds_test[var_key].max().item()}")
        else:
            logger.warning(f"Variable {var_key} not found in uxarray dataset!")
            logger.info(f"Available variables: {list(uxds_test.data_vars)}")
        logger.info("================================")
    else:
        logger.warning("uxds_test is None! Cannot plot native grid data without a grid file.")

    for region in regions:
        ds_test_region = _process_test_dataset(
            parameter, ds_test, ds_land_sea_mask, var_key, region
        )

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        print(f"Processing {var_key} for region {region}")

        # Use the simplified plotting function that works directly with uxarray datasets
        _save_plot(
            parameter,
            var_key,
            region,
            uxds_test=uxds_test,
            uxds_ref=None,
        )


def _run_diags_3d_model_only(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
):
    """Run a model-only diagnostics on a 3D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    """
    plev = parameter.plevs
    logger.info("Selected pressure level(s): {}".format(plev))

    ds_test_rg = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)

    for ilev in plev:
        z_axis_key = get_z_axis(ds_test_rg[var_key]).name
        ds_test_ilev = ds_test_rg.sel({z_axis_key: ilev})

        for region in regions:
            ds_test_region = _process_test_dataset(
                parameter, ds_test_ilev, ds_land_sea_mask, var_key, region
            )
            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(f"Processing {var_key} for region {region} at level {ilev}")

            # Use the simplified plotting function that works directly with uxarray datasets
            _save_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=None,
            )


def _check_var_dims(
    ds_test: xr.Dataset, ds_ref: xr.Dataset, var_key: str
) -> Tuple[bool, bool]:
    """Check if the variables have 3D dimensions and if their dimensions are different.

    Parameters
    ----------
    ds_test : xr.Dataset
        The test dataset.
    ds_ref : xr.Dataset
        The reference dataset.
    var_key : str
        The key of the variable.

    Returns
    -------
    Tuple[bool, bool]
        A tuple containing two boolean values:
        - is_vars_3d: True if both variables have 3D dimensions, False otherwise.
        - is_dims_diff: True if the dimensions of the two variables are different, False otherwise.
    """
    dv_test = ds_test[var_key]
    dv_ref = ds_ref[var_key]

    is_vars_3d = has_z_axis(dv_test) and has_z_axis(dv_ref)
    is_dims_diff = has_z_axis(dv_test) != has_z_axis(dv_ref)

    return is_vars_3d, is_dims_diff


def _run_diags_2d(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Run diagnostics on a 2D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    for region in regions:
        (
            ds_test_region,
            ds_test_region_regrid,
            ds_ref_region,
            ds_ref_region_regrid,
            ds_diff_region,
        ) = _subset_and_align_native_datasets(
            parameter,
            ds_test,
            ds_ref,
            ds_land_sea_mask,
            var_key,
            region,
            uxds_test,
            uxds_ref,
        )

        metrics_dict = _create_metrics_dict(
            var_key,
            ds_test_region,
            ds_test_region_regrid,
            ds_ref_region,
            ds_ref_region_regrid,
            ds_diff_region,
            uxds_test=uxds_test,
            uxds_ref=uxds_ref,
        )

        parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev=None)
        print(f"Processing {var_key} for region {region} (model vs obs)")

        # Use the simplified plotting function that works directly with uxarray datasets
        _save_plot(
            parameter,
            var_key,
            region,
            uxds_test=uxds_test,
            uxds_ref=uxds_ref,
        )


def _run_diags_3d(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    season: str,
    regions: List[str],
    var_key: str,
    ref_name: str,
    uxds_test: ux.dataset.UxDataset = None,
    uxds_ref: ux.dataset.UxDataset = None,
):
    """Run diagnostics on a 3D variable using native grid.

    This function gets the variable's metrics by region, then saves the
    metrics, metric plots, and data (optional, `CoreParameter.save_netcdf`).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the ref variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset, which is only used for masking if the region
        is "land" or "ocean".
    season : str
        The season.
    regions : List[str]
        The list of regions.
    var_key : str
        The key of the variable.
    ref_name : str
        The reference name.
    uxds_test : ux.dataset.UxDataset, optional
        The uxarray dataset containing the test native grid information.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid information.
    """
    plev = parameter.plevs
    logger.info("Selected pressure level(s): {}".format(plev))

    ds_test_rg = regrid_z_axis_to_plevs(ds_test, var_key, parameter.plevs)
    ds_ref_rg = regrid_z_axis_to_plevs(ds_ref, var_key, parameter.plevs)

    for ilev in plev:
        z_axis_key = get_z_axis(ds_test_rg[var_key]).name
        ds_test_ilev = ds_test_rg.sel({z_axis_key: ilev})
        ds_ref_ilev = ds_ref_rg.sel({z_axis_key: ilev})

        for region in regions:
            (
                ds_test_region,
                ds_test_region_regrid,
                ds_ref_region,
                ds_ref_region_regrid,
                ds_diff_region,
            ) = _subset_and_align_native_datasets(
                parameter,
                ds_test_ilev,
                ds_ref_ilev,
                ds_land_sea_mask,
                var_key,
                region,
                uxds_test,
                uxds_ref,
            )

            metrics_dict = _create_metrics_dict(
                var_key,
                ds_test_region,
                ds_test_region_regrid,
                ds_ref_region,
                ds_ref_region_regrid,
                ds_diff_region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
            )

            parameter._set_param_output_attrs(var_key, season, region, ref_name, ilev)
            print(f"Processing {var_key} for region {region} at level {ilev} (model vs obs)")

            # Use the simplified plotting function that works directly with uxarray datasets
            _save_plot(
                parameter,
                var_key,
                region,
                uxds_test=uxds_test,
                uxds_ref=uxds_ref,
            )


def _get_ref_climo_dataset(
    dataset: Dataset, var_key: str, season: ClimoFreq
) -> xr.Dataset | None:
    """Get the reference climatology dataset for the variable and season.

    If the reference climatatology does not exist or could not be found, it
    will be considered a model-only run and return `None`.

    Parameters
    ----------
    dataset : Dataset
        The dataset object.
    var_key : str
        The key of the variable.
    season : CLIMO_FREQ
        The climatology frequency.

    Returns
    -------
    xr.Dataset | None
        The reference climatology if it exists or None if it does not.
        None indicates a model-only run.

    Raises
    ------
    RuntimeError
        If `self.data_type` is not "ref".
    """
    if dataset.data_type == "ref":
        try:
            ds_ref = dataset.get_climo_dataset(var_key, season)
        except (RuntimeError, IOError):
            ds_ref = None

            logger.info("Cannot process reference data, analyzing test data only.")
    else:
        raise RuntimeError(
            "`Dataset._get_ref_dataset` only works with "
            f"`parameter.data_type == 'ref'`, not {dataset.data_type}."
        )

    return ds_ref


def _process_test_dataset(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
) -> xr.Dataset:
    """Process the test dataset for the given region.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The core parameter object.
    ds_test : xr.Dataset
        The test dataset.
    ds_land_sea_mask : xr.Dataset
        The dataset for land-sea mask.
    var_key : str
        The variable key.
    region : str
        The region to process.

    Returns:
    --------
    xr.Dataset
        The processed test dataset for the region.
    """
    ds_test_region = ds_test.copy()

    if "land" in region or "ocean" in region:
        ds_test_region = _apply_land_sea_mask(
            ds_test,
            ds_land_sea_mask,
            var_key,
            region,  # type: ignore
            parameter.regrid_tool,
            parameter.regrid_method,
        )

    if "global" not in region:
        ds_test_region = _subset_on_region(ds_test_region, var_key, region)

    return ds_test_region


def _subset_and_align_native_datasets(
    parameter: LatLonNativeParameter,
    ds_test: xr.Dataset,
    ds_ref: xr.Dataset,
    ds_land_sea_mask: xr.Dataset,
    var_key: str,
    region: str,
    uxds: ux.dataset.UxDataset = None,
):
    """Subset and align datasets using native grid data.

    For native grid data, this is an adapted version that preserves the native grid
    structure while still computing metrics and difference fields.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    ds_test : xr.Dataset
        The dataset containing the test variable.
    ds_ref : xr.Dataset
        The dataset containing the reference variable.
    ds_land_sea_mask : xr.Dataset
        The land sea mask dataset.
    var_key : str
        The variable key.
    region : str
        The region.
    uxds : ux.dataset.UxDataset, optional
        The uxarray dataset containing the native grid information.

    Returns
    -------
    Tuple[xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset, xr.Dataset]
        A tuple containing the processed datasets:
        - ds_test_region: The test dataset subsetted to the region.
        - ds_test_region_regrid: The test dataset regridded (same as ds_test_region for native).
        - ds_ref_region: The reference dataset subsetted to the region.
        - ds_ref_region_regrid: The reference dataset regridded to match the test data resolution.
        - ds_diff_region: The difference between the test and reference datasets.
    """
    # First, process the test and reference datasets for the given region
    ds_test_region = _process_test_dataset(
        parameter, ds_test, ds_land_sea_mask, var_key, region
    )

    ds_ref_region = _process_test_dataset(
        parameter, ds_ref, ds_land_sea_mask, var_key, region
    )

    # For native grid, we don't regrid the test data, we keep it on its native grid
    ds_test_region_regrid = ds_test_region

    # For model vs obs comparison, we still need to handle the reference data
    # For now, we'll use a simple approach - the reference data stays on its grid
    # but metrics calculations will handle the irregular grid
    ds_ref_region_regrid = ds_ref_region

    # For difference calculation, this is more complex with unstructured grids
    # In a first implementation, we'll compute a simple difference without regridding
    # In future versions, we might want to implement more sophisticated methods
    # Create a copy of the test dataset for the difference
    ds_diff_region = ds_test_region_regrid.copy()

    # For basic implementation, we'll just subtract the means for the difference plot
    test_mean = spatial_avg(ds_test_region_regrid, var_key)
    ref_mean = spatial_avg(ds_ref_region_regrid, var_key)
    diff = test_mean - ref_mean

    # Create a new variable in the difference dataset with a constant value
    # This is a simplification for the first implementation
    # In future versions, we may want to implement proper regridding between native and reference grids
    ds_diff_region[var_key].values = ds_diff_region[var_key].values * 0 + diff

    return (
        ds_test_region,
        ds_test_region_regrid,
        ds_ref_region,
        ds_ref_region_regrid,
        ds_diff_region,
    )


def _create_metrics_dict(
    var_key: str,
    ds_test: xr.Dataset,
    ds_test_regrid: xr.Dataset | None,
    ds_ref: xr.Dataset | None,
    ds_ref_regrid: xr.Dataset | None,
    ds_diff: xr.Dataset | None,
    uxds: ux.dataset.UxDataset = None,
) -> MetricsDict:
    """Calculate metrics using the variable in the datasets for native grid visualization.

    Metrics include min value, max value, spatial average (mean), standard
    deviation, correlation (pearson_r), and RMSE. The default value for
    optional metrics is None.

    Parameters
    ----------
    var_key : str
        The variable key.
    ds_test : xr.Dataset
        The test dataset.
    ds_test_regrid : xr.Dataset | None
        The regridded test Dataset. This arg will be None if a model only
        run is performed.
    ds_ref : xr.Dataset | None
        The optional reference dataset. This arg will be None if a model only
        run is performed.
    ds_ref_regrid : xr.Dataset | None
        The optional regridded reference dataset. This arg will be None if a
        model only run is performed.
    ds_diff : xr.Dataset | None
        The difference between ``ds_test_regrid`` and ``ds_ref_regrid`` if both
        exist. This arg will be None if a model only run is performed.
    uxds : ux.dataset.UxDataset, optional
        The uxarray dataset containing the native grid information.

    Returns
    -------
    Metrics
        A dictionary with the key being a string and the value being either
        a sub-dictionary (key is metric and value is float) or a string
        ("unit").
    """
    # Extract these variables for reuse.
    var_test = ds_test[var_key]

    # xarray.DataArray.min() and max() returns a `np.ndarray` with a single
    # int/float element. Using `.item()` returns that single element.
    metrics_dict: MetricsDict = {
        "test": {
            "min": var_test.min().item(),
            "max": var_test.max().item(),
            "mean": spatial_avg(ds_test, var_key),  # type: ignore
            "std": std(ds_test, var_key),
        },
        "unit": ds_test[var_key].attrs["units"],
    }
    metrics_dict = _set_default_metric_values(metrics_dict)

    # Additional native grid information
    if uxds is not None:
        try:
            grid_info = {
                "grid_type": getattr(uxds, "grid_type", "unknown"),
                "ne": getattr(uxds, "ne", ""),
                "npe": getattr(uxds, "npe", ""),
                "element_count": len(uxds.face) if hasattr(uxds, "face") else 0,
            }
            metrics_dict["grid_info"] = grid_info
        except Exception as e:
            logger.warning(f"Error adding grid info to metrics: {e}")

    if ds_ref is not None:
        var_ref = ds_ref[var_key]

        metrics_dict["ref"] = {
            "min": var_ref.min().item(),
            "max": var_ref.max().item(),
            "mean": spatial_avg(ds_ref, var_key),  # type: ignore
            "std": std(ds_ref, var_key),
        }

    if ds_test_regrid is not None and ds_ref_regrid is not None:
        var_test_regrid = ds_test_regrid[var_key]
        metrics_dict["test_regrid"] = {
            "min": var_test_regrid.min().item(),
            "max": var_test_regrid.max().item(),
            "mean": spatial_avg(ds_test_regrid, var_key),  # type: ignore
            "std": std(ds_test_regrid, var_key),
        }

        var_ref_regrid = ds_ref_regrid[var_key]
        metrics_dict["ref_regrid"] = {
            "min": var_ref_regrid.min().item(),
            "max": var_ref_regrid.max().item(),
            "mean": spatial_avg(ds_ref_regrid, var_key),  # type: ignore
            "std": std(ds_ref_regrid, var_key),
        }

        # For native grid, the correlation and RMSE calculations might need to be adapted
        # In the first implementation, we'll use global means for simplicity
        metrics_dict["misc"] = {
            "rmse": abs(
                metrics_dict["test_regrid"]["mean"] - metrics_dict["ref_regrid"]["mean"]
            ),
            "corr": 0.0,  # Placeholder - proper correlation would require regridding
        }

    # for model-only run
    if ds_test is not None and ds_ref_regrid is None:
        metrics_dict["test_regrid"] = metrics_dict["test"]

    if ds_diff is not None:
        var_diff = ds_diff[var_key]

        metrics_dict["diff"] = {
            "min": var_diff.min().item(),
            "max": var_diff.max().item(),
            "mean": spatial_avg(ds_diff, var_key),  # type: ignore
        }

    return metrics_dict


def _set_default_metric_values(metrics_dict: MetricsDict) -> MetricsDict:
    """Set the default values for the metrics in case they are not calculated.

    Parameters
    ----------
    metrics_dict : MetricsDict
        The dictionary containing the metrics.

    Returns
    -------
    MetricsDict
        The dictionary containing the metrics with default values.
    """
    var_keys = ["test_regrid", "ref", "ref_regrid", "diff"]
    metric_keys = ["min", "max", "mean", "std"]
    for var_key in var_keys:
        metrics_dict[var_key] = {
            metric_key: METRICS_DEFAULT_VALUE for metric_key in metric_keys
        }

    metrics_dict["misc"] = {
        "rmse": METRICS_DEFAULT_VALUE,
        "corr": METRICS_DEFAULT_VALUE,
    }

    return metrics_dict
