.. E3SM Diags documentation master file, created by
   sphinx-quickstart on Tue Sep  5 17:27:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. _index-label:

***************************
E3SM Diagnostics Package v3
***************************
Welcome to the E3SM Diagnostics Package documentation hub.

To change the documentation version, use the version selector in the bottom left-hand corner.
Please note, documentation for versions ``v2.5.0`` are not available in the version selector.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   input-data-requirement
   quickguides/index
   examples
   install
   config-run
   defining-parameters
   available-parameters
   colormaps
   add-new-diagnostics
   dev_guide/index
   contributing

Overview
--------

This diagnostics package is constructed to support the diagnostics
needs of DOE's `Energy Exascale Earth System Model (E3SM) project
<https://climatemodeling.science.energy.gov/projects/energy-exascale-earth-system-model>`__.
The ultimate goal of this work is to develop a comprehensive diagnostics package
that:

-  integrates the basic functionality of NCAR's AMWG diagnostics
   package;
-  utilizes most updated observational datasets, including remote
   sensing, reanalysis and in-situ datasets;
-  interfaces with diagnostics developed from different E3SM focus
   groups: atmosphere group, coupled simulation group, land group;
-  is flexible for user-specified diagnostics and configuration for
   use by other earth system models.

Current State
-------------

Algorithm and visualization codes for **latitude-longitude contour maps**,
**polar contour maps**, the accompanying **summarizing table** and
**Taylor diagram plots**,  **pressure-latitude zonal mean contour plots**,
**zonal mean line plots**, **pressure-longitude meridional mean contour plots**,
**area mean time series plots**, and **Cloud Top Height-Tau** joint histograms from
COSP cloud simulator output. Plots can be created for annual and seasonal climatologies,
and monthly mean time series. In additional to the core sets being released in v1,
**ENSO diags**, **QBO diags**, **Diurnal cycle phase plot**, **Streamflow evaluation**,
**ARM diags**, and **TC analysis**, **Mixed Phase Partition**, and **Wheeler-Kiladis Diagram** are implemented in v2 release. v3 introduces a
completely new back-end, replacing CDAT with Xarray and xCDAT. Due to the significant
scale of code changes, this has been incremented as a major release. The user-facing
API for running E3SM Diagnostics remains backward-compatible between v2 and v3.

The package also supports custom user diagnostics, by specifying
plot type, desired region (global, ocean, land, etc.),
pressure levels for variables with the vertical dimension.

For flexibility, the code structure cleanly separates data manipulation
(reading input files, processing data, etc) from plotting functions.
One back-end is available:

* `matplotlib <https://matplotlib.org>`_/ `cartopy <http://scitools.org.uk/cartopy>`_ (**mpl**)

Additional back-ends could be implemented if the need arose.

+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig01.png                    | .. figure:: _static/index/fig02.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig01.png                    |    :target: _static/index/fig02.png                  |
|                                                        |                                                      |
|    Latitude-longitude contour map                      |    Summary table                                     |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig03.png                    | .. figure:: _static/index/fig04.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig03.png                    |    :target: _static/index/fig04.png                  |
|                                                        |                                                      |
|    Taylor Diagram                                      |    Zonal mean line plot                              |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig05.png                    | .. figure:: _static/index/fig06.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig05.png                    |    :target: _static/index/fig06.png                  |
|                                                        |                                                      |
|    Pressure-latitude zonal mean contour                |    Polar contour map                                 |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig07.png                    | .. figure:: _static/index/fig08.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig07.png                    |    :target: _static/index/fig08.png                  |
|                                                        |                                                      |
|    Cloud Top Height-Tau joint histograms               |    Pressure-longitude meridional mean contour        |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig09.png                    | .. figure:: _static/index/fig10.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig09.png                    |    :target: _static/index/fig10.png                  |
|                                                        |                                                      |
|    ENSO diagnostics map                                |    ENSO diagnostics scatter plot                     |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig11.png                    | .. figure:: _static/index/fig12.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig11.png                    |    :target: _static/index/fig12.png                  |
|                                                        |                                                      |
|    QBO                                                 |    Area Mean time series                             |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig13.png                    | .. figure:: _static/index/fig14.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig13.png                    |    :target: _static/index/fig14.png                  |
|                                                        |                                                      |
|    Diurnal cycle phase maps                            |    Streamflow seasonality map                        |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig15.png                    | .. figure:: _static/index/fig16.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig15.png                    |    :target: _static/index/fig16.png                  |
|                                                        |                                                      |
|    Mean annual streamflow map                          |    Mean annual streamflow scatter plot               |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig17.png                    | .. figure:: _static/index/fig18.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig17.png                    |    :target: _static/index/fig18.png                  |
|                                                        |                                                      |
|    ARM diagnostics annual cycle plot                   |    ARM diagnostics diurnal cycle plot                |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig19.png                    | .. figure:: _static/index/fig20.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig19.png                    |    :target: _static/index/fig20.png                  |
|                                                        |                                                      |
|    ARM diagnostics monthly diurnal cycle of cloud plot |    ARM diagnostics convection onset statistics plot  |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig21.png                    | .. figure:: _static/index/fig22.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig21.png                    |    :target: _static/index/fig22.png                  |
|                                                        |                                                      |
|    Tropical Cyclone Track Density                      |    Annual Cycle Zonal Mean plot                      |
+--------------------------------------------------------+------------------------------------------------------+
| .. figure:: _static/index/fig23.png                    | .. figure:: _static/index/fig24.png                  |
|    :align: center                                      |    :align: center                                    |
|    :target: _static/index/fig23.png                    |    :target: _static/index/fig24.png                  |
|                                                        |                                                      |
|    Tropical Cyclone frequency per basin                |    Per-basin Tropical Cyclone frac seasonal cycle    |
+--------------------------------------------------------+------------------------------------------------------+

The above plots and more can be found
`here <https://portal.nersc.gov/cfs/e3sm/chengzhu/tutorial2024/e3sm_diags_extended/viewer/>`_.

Feature availability for each backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. |check| unicode:: U+2714  .. checkmark symbol
.. |ballot| unicode:: U+2718  .. ballot symbol

The table below summarizes current status of features supported by the backend.

+--------------------------------------------+---------+
| Plot set or Feature                        | mpl     |
+============================================+=========+
| Latitude-longitude contour maps            | |check| |
+--------------------------------------------+---------+
| Polar contour maps                         | |check| |
+--------------------------------------------+---------+
| Pressure-latitude zonal mean contour plots | |check| |
+--------------------------------------------+---------+
| Pressure-longitude meridional mean contour | |check| |
+--------------------------------------------+---------+
| Zonal mean line plots                      | |check| |
+--------------------------------------------+---------+
| Cloud Top Height-Tau joint histograms      | |check| |
+--------------------------------------------+---------+
| Area Mean time series plots                | |check| |
+--------------------------------------------+---------+
| Multi-processing                           | |check| |
+--------------------------------------------+---------+
