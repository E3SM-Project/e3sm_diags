from __future__ import annotations

from typing import TYPE_CHECKING, Optional

import cartopy.crs as ccrs
import cartopy.feature as cfeature
import matplotlib
import matplotlib.colors as mcolors
import numpy as np
import uxarray as ux

from e3sm_diags.derivations.default_regions_xr import REGION_SPECS
from e3sm_diags.logger import _setup_child_logger
from e3sm_diags.plot.utils import (
    DEFAULT_PANEL_CFG,
    _get_colormap,
    _save_plot,
)

if TYPE_CHECKING:
    from e3sm_diags.parameter.lat_lon_native_parameter import LatLonNativeParameter

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # isort:skip  # noqa: E402

logger = _setup_child_logger(__name__)


def plot(  # noqa: C901
    parameter: LatLonNativeParameter,
    var_key: str,
    region: str,
    uxds_test: ux.dataset.UxDataset,
    uxds_ref: Optional[ux.dataset.UxDataset] = None,
    uxds_diff: Optional[ux.dataset.UxDataset] = None,
):
    """Create visualization of data on native (unstructured) grids using uxarray.

    This function creates plots without regridding to a regular lat-lon grid.
    The layout matches the standard lat_lon_plot with 3 panels (test, ref, diff).

    Parameters
    ----------
    parameter : LatLonNativeParameter
        The parameter object.
    var_key : str
        The variable key.
    region : str
        The region name, used to determine map extents.
    uxds_test : ux.dataset.UxDataset
        The test native grid dataset.
    uxds_ref : ux.dataset.UxDataset, optional
        The reference native grid dataset.
    uxds_diff : ux.dataset.UxDataset, optional
        The difference native grid dataset.
    ilev : float, optional
        The pressure level to visualize for 3D variables.
    """
    logger.info(f"Creating native grid plot for {var_key}, region={region}")

    if uxds_test is None or var_key not in uxds_test:
        logger.error(
            f"Cannot plot native grid data. Either uxds_test is None or {var_key} not in dataset"
        )
        if uxds_test is not None:
            logger.error(f"Available variables: {list(uxds_test.data_vars)}")
        return

    has_reference = uxds_ref is not None and var_key in uxds_ref
    if has_reference:
        logger.info(
            f"Reference data available for {var_key}, implementing model vs model visualization"
        )

    # Set the viewer description to the "long_name" attr of the variable
    if "long_name" in uxds_test[var_key].attrs:
        parameter.viewer_descr[var_key] = uxds_test[var_key].attrs["long_name"]
    else:
        parameter.viewer_descr[var_key] = var_key

    # Create figure with standard layout
    fig = plt.figure(figsize=parameter.figsize, dpi=parameter.dpi)

    # Use the same title formatting as in lat_lon_plot (fontsize=18, y=0.96)
    fig.suptitle(parameter.main_title, x=0.5, y=0.96, fontsize=18)

    # Get region information for setting map extents
    region_specs = REGION_SPECS.get(region, None)

    # Set map bounds based on region
    if region_specs is not None:
        lat_bounds = region_specs.get("lat", (-90, 90))  # type: ignore
        lon_bounds = region_specs.get("lon", (0, 360))  # type: ignore
        is_global_domain = lat_bounds == (-90, 90) and lon_bounds == (0, 360)
    else:
        lat_bounds = (-90, 90)
        lon_bounds = (0, 360)
        is_global_domain = True

    # Get the cartopy projection based on region info.
    # --------------------------------------------------------------------------
    # Determine projection and extents based on region
    projection = ccrs.PlateCarree()
    if is_global_domain:
        projection = ccrs.PlateCarree(central_longitude=180)

    logger.info(f"Region: {region}, lat_bounds: {lat_bounds}, lon_bounds: {lon_bounds}")

    # Extract metrics directly from the uxarray dataset
    if uxds_test is not None and var_key in uxds_test:
        # ------------------------------------------------------------
        # FIXME: Metrics extraction for test, ref, and diff datasets are duplicated,
        # extract to a helper function.
        # ------------------------------------------------------------
        if is_global_domain:
            test_min = uxds_test[var_key].min().item()
            test_max = uxds_test[var_key].max().item()
            # For native grid, use simple mean as approximation
            test_mean = uxds_test[var_key].mean().item()
        else:
            test_subset = uxds_test[var_key].subset.bounding_box(
                lon_bounds,
                lat_bounds,
            )
            test_min = test_subset.min().item()
            test_max = test_subset.max().item()
            # For native grid, use simple mean as approximation
            test_mean = test_subset.mean().item()

        units = uxds_test[var_key].attrs.get("units", "")
        # ------------------------------------------------------------
    else:
        # This should not happen since we check earlier, but just in case
        logger.error(f"Missing test data for variable {var_key} in native grid dataset")

        return

    # Extract metrics for reference data if available
    ref_min = ref_max = ref_mean = diff_min = diff_max = diff_mean = None
    if has_reference and uxds_ref is not None:
        # ------------------------------------------------------------
        # FIXME: Metrics extraction for test, ref, and diff datasets are duplicated,
        # extract to a helper function.
        # ------------------------------------------------------------
        if is_global_domain:
            ref_min = uxds_ref[var_key].min().item()
            ref_max = uxds_ref[var_key].max().item()
            ref_mean = uxds_ref[var_key].mean().item()
        else:
            ref_subset = uxds_ref[var_key].subset.bounding_box(lon_bounds, lat_bounds)
            ref_min = ref_subset.min().item()
            ref_max = ref_subset.max().item()
            ref_mean = ref_subset.mean().item()
        # ------------------------------------------------------------

        ref_units = uxds_ref[var_key].attrs.get("units", "")

        # Check if units match between test and reference
        if units != ref_units:
            logger.warning(
                f"Units mismatch between test ({units}) and reference ({ref_units})"
            )

        # Calculate approximate metrics for difference if not already available
        if uxds_diff is not None and var_key in uxds_diff:
            # ------------------------------------------------------------
            # FIXME: Metrics extraction for test, ref, and diff datasets are duplicated,
            # extract to a helper function.
            # ------------------------------------------------------------
            if is_global_domain:
                diff_min = uxds_diff[var_key].min().item()
                diff_max = uxds_diff[var_key].max().item()
                diff_mean = uxds_diff[var_key].mean().item()
            else:
                diff_subset = uxds_diff[var_key].subset.bounding_box(
                    lon_bounds,
                    lat_bounds,
                )
                diff_min = diff_subset.min().item()
                diff_max = diff_subset.max().item()
                diff_mean = diff_subset.mean().item()
            # ------------------------------------------------------------

    # Create panels following the lat_lon_plot layout
    # Panel 1: Test data (always created)
    ax1 = fig.add_axes(DEFAULT_PANEL_CFG[0], projection=projection)

    # Use the standard title configuration from utils._configure_titles
    # Format: (years_text on left, main title in center, units on right)
    ax1.set_title(parameter.test_name_yrs, loc="left", fontdict={"fontsize": 9.5})
    ax1.set_title(parameter.test_title, fontdict={"fontsize": 11.5})
    ax1.set_title(units, loc="right", fontdict={"fontsize": 9.5})

    # Initialize ax2 and ax3 as None - they'll only be created when reference data exists
    ax2 = None
    ax3 = None

    # Only create panels 2 and 3 when reference data is available
    if has_reference:
        # Panel 2: Reference data
        ax2 = fig.add_axes(DEFAULT_PANEL_CFG[1], projection=projection)
        ax2.set_title(parameter.ref_name_yrs, loc="left", fontdict={"fontsize": 9.5})
        ax2.set_title(parameter.reference_title, fontdict={"fontsize": 11.5})
        ax2.set_title(units, loc="right", fontdict={"fontsize": 9.5})

        # Panel 3: Difference plot
        ax3 = fig.add_axes(DEFAULT_PANEL_CFG[2], projection=projection)
        ax3.set_title("", loc="left", fontdict={"fontsize": 9.5})
        ax3.set_title(parameter.diff_title, fontdict={"fontsize": 11.5})
        ax3.set_title(units, loc="right", fontdict={"fontsize": 9.5})

    # Configure map settings for all created panels
    panels = [ax for ax in [ax1, ax2, ax3] if ax is not None]
    _configure_map_panels(
        panels, region, region_specs, lat_bounds, lon_bounds, is_global_domain
    )

    # Set up periodic elements for all native grid plots
    periodic_elements = None
    if hasattr(parameter, "split_periodic_elements"):
        periodic_elements = "split" if parameter.split_periodic_elements else None

    # Create the test panel visualization
    _create_panel_visualization(
        uxds_test,
        var_key,
        ax1,
        fig,
        DEFAULT_PANEL_CFG[0],
        units,
        parameter.test_colormap,
        parameter.contour_levels,
        test_min,
        test_max,
        test_mean,
        False,
        periodic_elements,
        parameter.antialiased,
    )

    # Create reference panel visualization if available
    if has_reference and ax2 is not None:
        _create_panel_visualization(
            uxds_ref,
            var_key,
            ax2,
            fig,
            DEFAULT_PANEL_CFG[1],
            units,
            parameter.reference_colormap,
            parameter.contour_levels,
            ref_min,  # type: ignore
            ref_max,  # type: ignore
            ref_mean,  # type: ignore
            False,
            periodic_elements,
            parameter.antialiased,
        )

        # Create difference panel visualization if available
        if ax3 is not None:
            try:
                if uxds_diff is not None and var_key in uxds_diff:
                    _create_panel_visualization(
                        uxds_diff,
                        var_key,
                        ax3,
                        fig,
                        DEFAULT_PANEL_CFG[2],
                        units,
                        parameter.diff_colormap,
                        parameter.diff_levels
                        if hasattr(parameter, "diff_levels")
                        else None,
                        diff_min,  # type: ignore
                        diff_max,  # type: ignore
                        diff_mean,  # type: ignore
                        True,
                        periodic_elements,
                        parameter.antialiased,
                    )

                    # ------------------------------------------------------------
                    # FIXME: Re-use utils.add_rmse_corr_text here
                    # ------------------------------------------------------------
                    # Add RMSE and correlation text
                    rmse = (
                        abs(diff_mean) if diff_mean is not None else None
                    )  # Simplified RMSE
                    corr = 0.0  # Placeholder - proper correlation would require aligned data

                    if rmse is None:
                        logger.error("RMSE calculation failed: diff_mean is None")
                        return
                    elif rmse < 0.01:
                        rmse_fmt = "%.2e"
                    else:
                        rmse_fmt = "%.4f"

                    # Add the RMSE and CORR text using proper positions from utils._add_rmse_corr_text
                    # Position is set to (0.6635, -0.0105) for left text and (0.7635, -0.0105) for right text
                    fig.text(
                        DEFAULT_PANEL_CFG[2][0] + 0.6635,
                        DEFAULT_PANEL_CFG[2][1] - 0.0105,
                        "RMSE\nCORR",
                        ha="left",
                        fontdict={"fontsize": 9.5},
                    )

                    fig.text(
                        DEFAULT_PANEL_CFG[2][0] + 0.7635,
                        DEFAULT_PANEL_CFG[2][1] - 0.0105,
                        f"{rmse_fmt}\n%.4f" % (rmse, corr),
                        ha="right",
                        fontdict={"fontsize": 9.5},
                    )
                    # ------------------------------------------------------------

                else:
                    # If difference calculation failed, show a message
                    ax3.text(
                        0.5,
                        0.5,
                        "Could not calculate difference data for native grids",
                        transform=ax3.transAxes,
                        ha="center",
                        va="center",
                        fontsize=11,
                    )
            except Exception as e:
                # Fallback if there's an error in diff calculation or visualization
                logger.error(f"Error calculating or visualizing difference data: {e}")
                import traceback

                logger.error(traceback.format_exc())
                ax3.text(
                    0.5,
                    0.5,
                    f"Error in difference calculation:\n{str(e)}",
                    transform=ax3.transAxes,
                    ha="center",
                    va="center",
                    fontsize=11,
                )

    # Save the plot using the standard output path structure
    _save_plot(fig, parameter)
    plt.close(fig)


def _configure_map_panels(
    panels, region, region_specs, lat_bounds, lon_bounds, is_global_domain
):
    """Configure map settings (projection, extent, features) for all panels.

    FIXME: Refactor this function for readability.

    Parameters
    ----------
    panels : list
        List of matplotlib axes to configure
    region : str
        Region name
    region_specs : dict
        Region specifications from REGION_SPECS
    lat_bounds : tuple
        Latitude bounds (south, north)
    lon_bounds : tuple
        Longitude bounds (west, east)
    is_global_domain : bool
        Whether this is a global domain
    """
    # Determine X and Y ticks
    lat_south, lat_north = lat_bounds
    lon_west, lon_east = lon_bounds

    for ax in panels:
        # Handle global domain specially - don't use set_extent for global domain
        if is_global_domain:
            logger.info("Using global view")
            ax.set_global()
        else:
            try:
                # More robust longitude handling for map extents
                lon_west_orig, lon_east_orig = lon_west, lon_east

                # For regions that don't specify longitude (like 60S60N), use the full longitude range
                if region_specs and "lon" not in region_specs:
                    logger.info(
                        f"Region {region} only specifies latitude bounds, using full longitude range"
                    )
                    lon_west = 0
                    lon_east = 360

                # Now determine the best projection based on the region
                # For full longitude range or close to it, use central_longitude=180
                is_lon_full = lon_east - lon_west >= 350

                # Set up appropriate projection
                if is_lon_full:
                    logger.info(
                        "Using central longitude 180 for full/near-full longitude range"
                    )
                    projection = ccrs.PlateCarree(central_longitude=180)
                    ax.projection = projection
                    # For full longitude, use simplified extent setting
                    ax.set_extent(
                        [-180, 180, lat_south, lat_north], crs=ccrs.PlateCarree()
                    )
                else:
                    # For partial longitude ranges, we need to handle differently
                    # Normalize to [-180, 180] range for consistency with cartopy
                    if lon_west > 180:
                        lon_west -= 360
                    if lon_east > 180:
                        lon_east -= 360

                    # Handle cases where the region crosses the dateline
                    if lon_east < lon_west:
                        # This is a dateline-crossing region (e.g., Pacific)
                        logger.info(
                            f"Detected dateline crossing region: lon=[{lon_west}, {lon_east}]"
                        )

                        # For dateline-crossing regions, adjust the central longitude of projection
                        center_lon = (lon_west + lon_east + 360) / 2.0
                        if center_lon > 180:
                            center_lon -= 360

                        logger.info(f"Using central longitude: {center_lon}")

                        # Create a new projection with the adjusted central longitude
                        ax.projection = ccrs.PlateCarree(central_longitude=center_lon)

                        # When using a central_longitude, we need to transform our coordinates
                        # Adjust longitudes for the new central longitude
                        if lon_west < 0:
                            lon_west += 360
                        if lon_east < 0:
                            lon_east += 360

                        logger.info(
                            f"Transformed coordinates for central_longitude={center_lon}: lon=[{lon_west}, {lon_east}]"
                        )

                    # Make sure longitudes are properly ordered
                    if lon_east < lon_west:
                        logger.warning(
                            "East longitude still less than west after transforms - swapping values"
                        )
                        lon_west, lon_east = lon_east, lon_west

                    logger.info(
                        f"Final map extent: lon=[{lon_west}, {lon_east}], lat=[{lat_south}, {lat_north}]"
                    )

                    # Set the extent using the adjusted longitude values
                    ax.set_extent(
                        [lon_west, lon_east, lat_south, lat_north],
                        crs=ccrs.PlateCarree(),
                    )

            except Exception as e:
                # Comprehensive error handling
                logger.error(f"Error setting map extent: {e}")
                logger.error(f"Original lon bounds: [{lon_west_orig}, {lon_east_orig}]")
                logger.error(f"Transformed lon bounds: [{lon_west}, {lon_east}]")
                import traceback

                logger.error(traceback.format_exc())

                # Fallback to global view
                logger.info("Falling back to global view due to extent error")
                ax.set_global()

        # Add map features
        ax.coastlines(linewidth=0.5)
        ax.add_feature(cfeature.BORDERS, linewidth=0.3)

        # Configure gridlines and labels
        gl = ax.gridlines(
            crs=ccrs.PlateCarree(),
            draw_labels=True,
            linewidth=0.5,
            color="gray",
            alpha=0.5,
            linestyle="--",
        )
        gl.top_labels = False
        gl.right_labels = False

        # Set aspect ratio only for non-global views
        if not is_global_domain:
            # Ensure we have a valid aspect ratio by clamping to reasonable values
            width = lon_east - lon_west
            height = lat_north - lat_south
            if width <= 0:
                width += 360  # Handle wraparound cases
            aspect_ratio = width / (
                2 * max(height, 1)
            )  # Avoid division by zero or negative values
            logger.info(f"Setting aspect ratio: {aspect_ratio}")
            ax.set_aspect(aspect_ratio)


def _create_panel_visualization(
    dataset: ux.UxDataset,
    var_key: str,
    ax: matplotlib.axes.Axes,
    fig: matplotlib.figure.Figure,
    panel_cfg: tuple[float, float, float, float],
    units: str,
    colormap_name: str,
    contour_levels: list[float] | None,
    min_value: float,
    max_value: float,
    mean_value: float,
    is_diff: bool = False,
    periodic_elements: str | None = None,
    antialiased: bool = True,
) -> matplotlib.collections.PolyCollection:
    """Create a panel visualization with PolyCollection for native grid data.

    Parameters
    ----------
    dataset : ux.UxDataset
        The uxarray dataset
    var_key : str
        The variable key
    ax : matplotlib.axes.Axes
        The axis to draw on
    fig : matplotlib.figure.Figure
        The figure for adding colorbars and text
    panel_cfg : tuple
        The panel configuration (x, y, width, height)
    units : str
        The units string
    colormap_name : str
        The name of the colormap
    contour_levels : list or None
        List of contour levels or None
    min_value, max_value, mean_value : float
        The min, max, and mean values
    is_diff : bool
        Whether this is a difference plot
    periodic_elements : str or None
        How to handle periodic elements in the PolyCollection
    antialiased : bool
        Whether to antialias the PolyCollection

    Returns
    -------
    pc : matplotlib.collections.PolyCollection
        The created PolyCollection
    """
    # Get the data array and handle time dimension if present
    var_data = dataset[var_key]

    # Check if time dimension exists and has more than one point
    if "time" in var_data.dims and var_data.sizes["time"] > 1:
        logger.warning(
            f"Variable {var_key} has multiple time points. Using first time point only."
        )
        # Select first time point
        var_data = var_data.isel(time=0)

    # Squeeze to remove any remaining singleton dimensions
    var_data = var_data.squeeze()

    # Log shape information for debugging
    logger.info(f"Variable {var_key} shape: {var_data.shape}")
    logger.info(f"Variable {var_key} dims: {var_data.dims}")

    # ------------------------------------------------------------
    # FIXME: Re-use utils._get_c_levels_and_norm here
    # ------------------------------------------------------------
    # Get colormap
    cmap = _get_colormap(colormap_name)

    # Create normalization
    if contour_levels and len(contour_levels) > 0:
        # Use the contour levels to create a BoundaryNorm
        levels = contour_levels
        # Create boundary norm with extended boundaries for values beyond the levels
        boundaries = [-1.0e8] + levels + [1.0e8]
        norm = mcolors.BoundaryNorm(boundaries=boundaries, ncolors=256)
    else:
        # For difference plots, use symmetric normalization
        if is_diff:
            max_abs = max(abs(min_value), abs(max_value))
            vmin, vmax = -max_abs, max_abs
        else:
            vmin, vmax = min_value, max_value

        # Add buffer for constant values
        if vmin == vmax:
            buffer = max(0.1, abs(vmin * 0.1))
            vmin -= buffer
            vmax += buffer
            logger.warning(f"Data has constant value, adding buffer: [{vmin}, {vmax}]")

        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)

    # Create the PolyCollection
    pc = var_data.to_polycollection(periodic_elements=periodic_elements)

    # Set visualization properties
    pc.set_cmap(cmap)
    pc.set_norm(norm)
    pc.set_antialiased(antialiased)

    # Add to panel
    ax.add_collection(pc)

    # ------------------------------------------------------------
    # FIXME: Re-use utils.add_colorbar here
    # ------------------------------------------------------------
    # Add colorbar
    cbax_rect = (
        panel_cfg[0] + 0.6635,  # Position relative to the panel
        panel_cfg[1] + 0.0215,
        0.0326,  # Width
        0.1792,  # Height
    )
    cbax = fig.add_axes(cbax_rect)
    cbar = fig.colorbar(pc, cax=cbax, extend="both")

    # Configure colorbar ticks
    if contour_levels and len(contour_levels) > 0:
        cbar.set_ticks(contour_levels)

        # Format tick labels
        maxval = np.amax(np.absolute(contour_levels))
        if maxval < 0.01:
            fmt, pad = "%.1e", 35
        elif maxval < 0.2:
            fmt, pad = "%5.3f", 28
        elif maxval < 10.0:
            fmt, pad = "%5.2f", 25
        elif maxval < 100.0:
            fmt, pad = "%5.1f", 25
        elif maxval > 9999.0:
            fmt, pad = "%.0f", 40
        else:
            fmt, pad = "%6.1f", 30

        labels = [fmt % level for level in contour_levels]
        cbar.ax.set_yticklabels(labels, ha="right")
        cbar.ax.tick_params(labelsize=9.0, pad=pad, length=0)
    else:
        cbar.ax.tick_params(labelsize=9.0, length=0)

    # Add units label
    cbar.set_label(units, fontsize=9.5)

    # Add metrics text (min, max, mean)
    metrics = (max_value, mean_value, min_value)

    # Format metrics based on size
    fmt_m = []
    for i in range(3):
        val = metrics[i]
        if abs(val) < 0.01 and val != 0:
            fs = "e"
        elif abs(val) > 100000.0:
            fs = "e"
        else:
            fs = "2f"
        fmt_m.append(fs)

    fmt_metrics = f"%.{fmt_m[0]}\n%.{fmt_m[1]}\n%.{fmt_m[2]}"

    # ------------------------------------------------------------
    # FIXME: Re-use utils.add_min_mean_max_text here
    # ------------------------------------------------------------
    # Add "Max, Mean, Min" labels
    fig.text(
        panel_cfg[0] + 0.6635,
        panel_cfg[1] + 0.2107,
        "Max\nMean\nMin",
        ha="left",
        fontdict={"fontsize": 9.5},
    )

    # Add the actual metric values
    fig.text(
        panel_cfg[0] + 0.7635,
        panel_cfg[1] + 0.2107,
        fmt_metrics % metrics,
        ha="right",
        fontdict={"fontsize": 9.5},
    )

    return pc
