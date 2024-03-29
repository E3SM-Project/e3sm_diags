"""Manages plotting, provides a single interface
for different plots with different backends."""
from __future__ import absolute_import, print_function

import importlib
import os
import sys
import traceback

import numpy
from matplotlib.colors import LinearSegmentedColormap

import e3sm_diags
from e3sm_diags.logger import custom_logger

logger = custom_logger(__name__)


def _get_plot_fcn(backend, set_name):
    """Get the actual plot() function based on the backend and set_name."""
    try:
        if backend in ["matplotlib", "mpl"]:
            backend = "cartopy"

        mod_str = "e3sm_diags.plot.{}.{}_plot".format(backend, set_name)
        module = importlib.import_module(mod_str)
        return module.plot

    except ModuleNotFoundError:
        logger.error(
            "Plotting for set {} with {} is not supported".format(set_name, backend)
        )
        traceback.print_exc()


def plot(set_name, ref, test, diff, metrics_dict, parameter):
    """Based on set_name and parameter.backend, call the correct plotting function.

    #TODO: Make metrics_dict a kwarg and update the other plot() functions
    """
    if hasattr(parameter, "plot"):
        parameter.plot(ref, test, diff, metrics_dict, parameter)
    else:
        if parameter.backend not in ["cartopy", "mpl", "matplotlib"]:
            raise RuntimeError('Invalid backend, use "matplotlib"/"mpl"/"cartopy"')

        plot_fcn = _get_plot_fcn(parameter.backend, set_name)
        if plot_fcn:
            try:
                plot_fcn(ref, test, diff, metrics_dict, parameter)
            except Exception as e:
                logger.exception(
                    "Error while plotting {} with backend {}".format(
                        set_name, parameter.backend
                    ),
                    exc_info=True,
                )
                traceback.print_exc()
                if parameter.debug:
                    sys.exit()


def get_colormap(colormap, parameters):
    """Get the colormap (string or mpl colormap obj), which can be
    loaded from a local file in the cwd, installed file, or a predefined mpl one."""
    colormap = str(colormap)  # unicode don't seem to work well with string.endswith()
    if not colormap.endswith(".rgb"):  # predefined vcs/mpl colormap
        return colormap

    installed_colormap = os.path.join(e3sm_diags.INSTALL_PATH, "colormaps", colormap)

    if os.path.exists(colormap):
        # colormap is an .rgb in the current directory
        pass
    elif not os.path.exists(colormap) and os.path.exists(installed_colormap):
        # use the colormap from /plot/colormaps
        colormap = installed_colormap
    elif not os.path.exists(colormap) and not os.path.exists(installed_colormap):
        pth = os.path.join(e3sm_diags.INSTALL_PATH, "colormaps")
        msg = "File {} isn't in the current working directory or installed in {}"
        raise IOError(msg.format(colormap, pth))

    rgb_arr = numpy.loadtxt(colormap)
    rgb_arr = rgb_arr / 255.0

    if parameters.backend in ["cartopy", "mpl", "matplotlib"]:
        cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)
        return cmap

    else:
        raise RuntimeError("Invalid backend: {}".format(parameters.backend))
