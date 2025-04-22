from __future__ import annotations

from typing import TYPE_CHECKING, List, Tuple

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import uxarray as ux
import xarray as xr

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    _add_colorbar,
    _add_min_mean_max_text,
    _add_rmse_corr_text,
    _configure_titles,
    _configure_x_and_y_axes,
    _get_c_levels_and_norm,
    _get_x_ticks,
    _get_y_ticks,
    _save_plot,
)

if TYPE_CHECKING:
    from e3sm_diags.driver.lat_lon_native_driver import MetricsDict
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)


def plot(
    parameter: LatLonNativeParameter,
    metrics_dict: MetricsDict,
    uxds_test: ux.dataset.UxDataset,
    uxds_ref: ux.dataset.UxDataset = None,
    uxds_diff: ux.dataset.UxDataset = None,
):
    """Plot the variable's metrics generated for the lat_lon_native set.

    This function creates visualizations directly from uxarray datasets
    without relying on regridded xarray DataArrays.

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object containing plot configurations.
    metrics_dict : MetricsDict
        The metrics calculated for all datasets.
    uxds_test : ux.dataset.UxDataset
        The uxarray dataset containing the test native grid data.
    uxds_ref : ux.dataset.UxDataset, optional
        The uxarray dataset containing the reference native grid data.
    uxds_diff : ux.dataset.UxDataset, optional
        The uxarray dataset containing the difference native grid data.
    """
    # Get the variable key from parameter
    var_key = parameter.var_id

    # Check that the test data contains the variable
    if var_key not in uxds_test:
        logger.error(f"Variable {var_key} not found in test dataset!")
        logger.info(f"Available variables: {list(uxds_test.data_vars)}")
        return

    # Create the figure for the plots
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Get the variable units
    units = metrics_dict["unit"]
    logger.info(f"Using units: {units} for variable {var_key}")

    # Extract metrics for test data
    min1 = metrics_dict["test"]["min"]  # type: ignore
    mean1 = metrics_dict["test"]["mean"]  # type: ignore
    max1 = metrics_dict["test"]["max"]  # type: ignore

    # Log important information about the test dataset
    logger.info(f"Plotting test data for variable {var_key}")
    logger.info(f"Test data metrics - min: {min1}, mean: {mean1}, max: {max1}")

    # Check for potential issues with the test data
    if min1 == max1:
        logger.warning(f"⚠️ Test data has identical min and max values: {min1}")

    # Create the test panel using native grid data
    _add_native_colormap(
        0,  # First panel
        var_key,  # Use the variable key instead of DataArray
        fig,
        parameter,
        parameter.test_colormap,
        parameter.contour_levels,
        title=(parameter.test_name_yrs, parameter.test_title, units),  # type: ignore
        metrics=(max1, mean1, min1),  # type: ignore
        uxds=uxds_test,  # Pass the test native grid dataset
    )

    # Add reference and difference panels if available
    if uxds_ref is not None and var_key in uxds_ref:
        # Extract metrics for reference data
        min2 = metrics_dict["ref"]["min"]  # type: ignore
        mean2 = metrics_dict["ref"]["mean"]  # type: ignore
        max2 = metrics_dict["ref"]["max"]  # type: ignore

        logger.info(f"Plotting reference data for variable {var_key}")
        logger.info(f"Reference data metrics - min: {min2}, mean: {mean2}, max: {max2}")

        # Create the reference panel using native grid data if available
        _add_native_colormap(
            1,  # Second panel
            var_key,  # Use the variable key instead of DataArray
            fig,
            parameter,
            parameter.reference_colormap,
            parameter.contour_levels,
            title=(parameter.ref_name_yrs, parameter.reference_title, units),  # type: ignore
            metrics=(min2, mean2, max2),  # type: ignore
            uxds=uxds_ref,  # Pass the reference native grid dataset
        )

        # Check if we have difference data
        if uxds_diff is not None and var_key in uxds_diff:
            # Extract metrics for difference data
            min3 = metrics_dict["diff"]["min"]  # type: ignore
            mean3 = metrics_dict["diff"]["mean"]  # type: ignore
            max3 = metrics_dict["diff"]["max"]  # type: ignore
            r = metrics_dict["misc"]["rmse"]  # type: ignore
            c = metrics_dict["misc"]["corr"]  # type: ignore

            logger.info(f"Plotting difference data for variable {var_key}")
            logger.info(f"Difference metrics - min: {min3}, mean: {mean3}, max: {max3}, rmse: {r}, corr: {c}")

            # Create the difference panel using native grid data
            _add_native_colormap(
                2,  # Third panel
                var_key,  # Use the variable key instead of DataArray
                fig,
                parameter,
                parameter.diff_colormap,
                parameter.diff_levels,
                title=(None, parameter.diff_title, units),  # type: ignore
                metrics=(max3, mean3, min3, r, c),  # type: ignore
                uxds=uxds_diff,  # Pass the difference native grid dataset
            )
        else:
            # If no difference dataset is available, create an empty panel with a message
            logger.info("No difference dataset available, creating empty difference panel")
            _create_empty_panel(
                2,  # Third panel
                fig,
                parameter,
                "No difference data available",
                units,
            )
    else:
        # If no reference dataset is available, create empty panels
        logger.info("No reference dataset available, creating empty panels")
        _create_empty_panel(
            1,  # Second panel
            fig,
            parameter,
            "No reference data available",
            units,
        )
        _create_empty_panel(
            2,  # Third panel
            fig,
            parameter,
            "No difference data available",
            units,
        )

    _save_plot(fig, parameter)

    plt.close()


def _add_native_colormap(
    subplot_num: int,
    var_key: str,
    fig: plt.Figure,
    parameter: LatLonNativeParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[str | None, str, str],
    metrics: Tuple[float, ...],
    uxds: ux.dataset.UxDataset = None,
):
    """Add a native grid colormap containing the variable data and metrics to the figure.

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    var : xr.DataArray
        The variable to plot.
    fig : plt.Figure
        The figure object to add the subplot to.
    parameter : LatLonNativeParameter
        The parameter object containing plot configurations.
    color_map : str
        The colormap styling to use (e.g., "cet_rainbow.rgb").
    contour_levels : List[float]
        The map contour levels.
    title : Tuple[str | None, str, str]
        A tuple of strings to form the title of the colormap, in the format
        (<optional> years, title, units).
    metrics : Tuple[float, ...]
        A tuple of metrics for this subplot.
    uxds : ux.dataset.UxDataset, optional
        The uxarray dataset containing the native grid information.
    """
    # Configure contour levels and boundary norm.
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Get region info and X and Y plot ticks.
    region_key = parameter.regions[0]
    region_specs = REGION_SPECS[region_key]

    # Get the region's domain slices for latitude and longitude if set, or
    # use the default value. If both are not set, then the region type is
    # considered "global".
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (0, 360))  # type: ignore

    # Boolean flags for configuring plots.
    is_global_domain = lat_slice == (-90, 90) and lon_slice == (0, 360)
    is_lon_full = lon_slice == (0, 360)

    # Determine X and Y ticks using longitude and latitude domains respectively.
    lon_west, lon_east = lon_slice
    x_ticks = _get_x_ticks(lon_west, lon_east, is_global_domain, is_lon_full)

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north)

    # Get the cartopy projection based on region info.
    projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        projection = ccrs.PlateCarree(central_longitude=180)

    # Get the figure Axes object using the projection above.
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=projection)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=projection)

    # For native grid plotting, use the uxarray PolyCollection approach
    if uxds is not None:
        try:
            # Use var_key directly instead of var.name
            if var_key in uxds:
                ux_var = uxds[var_key].squeeze()
                # Create the PolyCollection for plotting
                periodic_elements = (
                    "split" if parameter.split_periodic_elements else None
                )
                pc = ux_var.to_polycollection(periodic_elements=periodic_elements)

                # Set styling options
                pc.set_antialiased(parameter.antialiased)
                pc.set_cmap(color_map)
                pc.set_norm(norm)

                # Add edge styling if specified
                if parameter.edge_color is not None:
                    pc.set_edgecolor(parameter.edge_color)
                    pc.set_linewidth(parameter.edge_width)
                else:
                    pc.set_edgecolor("none")

                # Add the collection to the plot
                logger.info(pc.get_fc())
                ax.add_collection(pc)
                logger.info("Successfully plotted native grid data")
            else:
                logger.error(f"Variable {var_key} not found in UxDataset")
                # Add error text to the plot
                ax.text(0.5, 0.5, f"Variable {var_key} not found in dataset",
                       transform=ax.transAxes, ha='center', va='center',
                       fontsize=10, color='red')

                ## Use a dummy mappable for the colorbar
                #pc = plt.cm.ScalarMappable(cmap=color_map, norm=norm)
                #pc.set_array(np.linspace(0, 1, 10))
        except Exception as e:
            logger.error(f"Error in native grid plotting: {e}")
            # Add error text to the plot
            ax.text(0.5, 0.5, f"Error plotting native grid data:\n{str(e)}",
                   transform=ax.transAxes, ha='center', va='center',
                   fontsize=10, color='red')

            # pc will be created later for the colorbar
    else:
        logger.error("No native grid data (uxds) provided")
        # Add error text to the plot
        ax.text(0.5, 0.5, "No native grid dataset provided",
               transform=ax.transAxes, ha='center', va='center',
               fontsize=10, color='red')

        # pc will be created later for the colorbar

    # Configure the aspect ratio and coast lines.
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    if not is_global_domain and "RRM" in region_key:
        ax.coastlines(resolution="50m", color="black", linewidth=1)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Configure the titles, x and y axes, and colorbar.
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )

    # ADD THIS TEXT TO PLOT FOR DEBUGGING
    ax.text(0.5, 0.5, "COLORBAR DISABLED FOR TESTING",
           transform=ax.transAxes, ha='center', va='center',
           fontsize=14, color='red', weight='bold')

    logger.info("ALL COLORBAR CODE COMPLETELY REMOVED FOR TESTING")
    # No colorbar code at all

    # Add metrics text to the figure.
    _add_min_mean_max_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    # Add grid information
    if uxds is not None:
        _add_native_grid_info(fig, subplot_num, uxds, DEFAULT_PANEL_CFG)
    else:
        # Fallback for regular grid info
        _add_standard_grid_info(fig, subplot_num, region_key, var, DEFAULT_PANEL_CFG)


def _plot_fallback(ax, var, color_map, projection, norm, c_levels):
    """Fallback to regular plotting when native grid plotting fails."""
    import xcdat as xc

    # Make the longitude cyclic if needed
    from e3sm_diags.plot.utils import _add_contour_plot, _make_lon_cyclic

    var = _make_lon_cyclic(var)

    # Get the coordinates
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    _add_contour_plot(ax, var, lon, lat, color_map, ccrs.PlateCarree(), norm, c_levels)


def _add_standard_colormap(
    subplot_num: int,
    var: xr.DataArray,
    fig: plt.Figure,
    parameter: LatLonNativeParameter,
    color_map: str,
    contour_levels: List[float],
    title: Tuple[str | None, str, str],
    metrics: Tuple[float, ...],
):
    """Add a standard colormap for non-native grid data.

    This is used for reference and difference plots that are not on the native grid.
    Adapted from lat_lon_plot._add_colormap.
    """
    import xcdat as xc

    from e3sm_diags.plot.utils import _add_contour_plot, _make_lon_cyclic

    var = _make_lon_cyclic(var)
    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    var = var.squeeze()

    # Configure contour levels and boundary norm.
    c_levels, norm = _get_c_levels_and_norm(contour_levels)

    # Get region info and X and Y plot ticks.
    region_key = parameter.regions[0]
    region_specs = REGION_SPECS[region_key]

    # Get the region's domain slices for latitude and longitude if set, or
    # use the default value. If both are not set, then the region type is
    # considered "global".
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (0, 360))  # type: ignore

    # Boolean flags for configuring plots.
    is_global_domain = lat_slice == (-90, 90) and lon_slice == (0, 360)
    is_lon_full = lon_slice == (0, 360)

    # Determine X and Y ticks using longitude and latitude domains respectively.
    lon_west, lon_east = lon_slice
    x_ticks = _get_x_ticks(lon_west, lon_east, is_global_domain, is_lon_full)

    lat_south, lat_north = lat_slice
    y_ticks = _get_y_ticks(lat_south, lat_north)

    # Get the cartopy projection based on region info.
    projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        projection = ccrs.PlateCarree(central_longitude=180)

    # Get the figure Axes object using the projection above.
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=projection)
    ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=projection)
    contour_plot = _add_contour_plot(
        ax, var, lon, lat, color_map, ccrs.PlateCarree(), norm, c_levels
    )

    # Configure the aspect ratio and coast lines.
    ax.set_aspect((lon_east - lon_west) / (2 * (lat_north - lat_south)))
    ax.coastlines(lw=0.3)

    if not is_global_domain and "RRM" in region_key:
        ax.coastlines(resolution="50m", color="black", linewidth=1)
        state_borders = cfeature.NaturalEarthFeature(
            category="cultural",
            name="admin_1_states_provinces_lakes",
            scale="50m",
            facecolor="none",
        )
        ax.add_feature(state_borders, edgecolor="black")

    # Configure the titles, x and y axes, and colorbar.
    _configure_titles(ax, title)
    _configure_x_and_y_axes(
        ax, x_ticks, y_ticks, ccrs.PlateCarree(), parameter.current_set
    )
    _add_colorbar(fig, subplot_num, DEFAULT_PANEL_CFG, contour_plot, c_levels)

    # Add metrics text to the figure.
    _add_min_mean_max_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    if len(metrics) == 5:
        _add_rmse_corr_text(fig, subplot_num, DEFAULT_PANEL_CFG, metrics)

    _add_standard_grid_info(fig, subplot_num, region_key, var, DEFAULT_PANEL_CFG)


def _add_native_grid_info(fig, subplot_num, uxds, panel_configs):
    """Add native grid information to the plot."""
    # Try to extract grid information from the uxarray dataset
    try:
        grid_type = getattr(uxds, "grid_type", "native")
        ne_val = getattr(uxds, "ne", "")
        if ne_val:
            grid_info = f"{grid_type} (ne={ne_val})"
        else:
            grid_info = f"{grid_type}"

        element_count = len(uxds.face) if hasattr(uxds, "face") else 0
        if element_count > 0:
            grid_info += f", {element_count} elements"
    except Exception:
        grid_info = "Native grid"

    # Get the panel config for the current subplot
    panel = panel_configs[subplot_num]

    # Create text position at the bottom-right of the subplot
    text_x = panel[0] + panel[2] - 0.01
    text_y = panel[1] + 0.01

    # Add the grid info text
    fig.text(
        text_x,
        text_y,
        grid_info,
        ha="right",
        va="bottom",
        color="black",
        fontsize=8,
    )


def _create_empty_panel(
    subplot_num: int,
    fig: plt.Figure,
    parameter: LatLonNativeParameter,
    message: str,
    units: str,
):
    """Create an empty panel with a message.

    Parameters
    ----------
    subplot_num : int
        The subplot number.
    fig : plt.Figure
        The figure object.
    parameter : LatLonNativeParameter
        The parameter object.
    message : str
        The message to display in the panel.
    units : str
        The units string (not used but included for compatibility).
    """
    logger.info(f"Creating empty panel {subplot_num} with message: {message}")

    # Get region info for consistent plotting
    region_key = parameter.regions[0]
    region_specs = REGION_SPECS[region_key]
    lat_slice = region_specs.get("lat", (-90, 90))  # type: ignore
    lon_slice = region_specs.get("lon", (0, 360))  # type: ignore
    is_global_domain = lat_slice == (-90, 90) and lon_slice == (0, 360)
    is_lon_full = lon_slice == (0, 360)

    # Create the axes with appropriate projection
    projection = ccrs.PlateCarree()
    if is_global_domain or is_lon_full:
        projection = ccrs.PlateCarree(central_longitude=180)

    # Add the axes to the figure
    ax = fig.add_axes(DEFAULT_PANEL_CFG[subplot_num], projection=projection)

    # Set extent based on region
    lon_west, lon_east = lon_slice
    lat_south, lat_north = lat_slice

    if is_global_domain:
        ax.set_global()
    else:
        try:
            ax.set_extent([lon_west, lon_east, lat_south, lat_north], crs=ccrs.PlateCarree())
        except Exception as e:
            logger.error(f"Error setting extent: {e}")
            ax.set_global()

    # Add basic map features
    ax.coastlines(linewidth=0.5)
    ax.add_feature(cfeature.BORDERS, linewidth=0.3)

    # Add the message text
    ax.text(0.5, 0.5, message,
           transform=ax.transAxes, ha='center', va='center',
           fontsize=12, style='italic')

    # No colorbar for empty panels
    logger.info("Empty panel created without colorbar")


def _add_standard_grid_info(fig, subplot_num, region_key, var, panel_configs):
    """Add standard (regular grid) information to the plot."""
    # Adapted from _add_grid_res_info in utils.py
    import xcdat as xc

    lat = xc.get_dim_coords(var, axis="Y")
    lon = xc.get_dim_coords(var, axis="X")

    # Get the number of lat/lon points
    try:
        num_lat = len(lat)
        num_lon = len(lon)

        # Calculate resolution
        try:
            lat_res = round(abs(lat[1] - lat[0]), 2)
            lon_res = round(abs(lon[1] - lon[0]), 2)
            res_text = f"{lat_res}°×{lon_res}°"
        except (IndexError, TypeError):
            res_text = ""

        grid_info = f"{region_key}: {num_lat}×{num_lon}"
        if res_text:
            grid_info += f" ({res_text})"
    except (TypeError, IndexError):
        grid_info = region_key

    # Get the panel config for the current subplot
    panel = panel_configs[subplot_num]

    # Create text position at the bottom-right of the subplot
    text_x = panel[0] + panel[2] - 0.01
    text_y = panel[1] + 0.01

    # Add the grid info text
    fig.text(
        text_x,
        text_y,
        grid_info,
        ha="right",
        va="bottom",
        color="black",
        fontsize=8,
    )
