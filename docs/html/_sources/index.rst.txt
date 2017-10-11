.. ACME Diags documentation master file, created by
   sphinx-quickstart on Tue Sep  5 17:27:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

************************
E3SM Diagnostics Package
************************

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   quick-guide-aims4
   quick-guide-edison
   install-config-run
   available-parameters
   colormaps
   add-new-diagnostics
   contributing


Overview
--------

This diagnostics package is constructed to support the diagnostics
needs of DOE's `Energy Exascale Earth System Model (E3SM) project 
<https://climatemodeling.science.energy.gov/projects/energy-exascale-earth-system-model>`__,
formerly known as Accelerated Climate Modeling for Energy (ACME).
The ultimate goal of this work is to develop a comprehensive diagnostics package
that:

-  fully integrates the functionality of NCAR's AMWG diagnostics
   package;
-  utilizes most updated observational datasets, including remote
   sensing, reanalysis and in-situ datasets;
-  interfaces with diagnostics developed from different E3SM focus
   groups: atmosphere group, coupled simulation group, land group;
-  interacts effectively with the PCMDI's metrics package and the ARM
   diagnostics package through a unifying framework: `Community
   Diagnostics Package (CDP) <https://github.com/UV-CDAT/cdp>`_.
-  is flexible for user-specified diagnostics and configuration for
   use by other climate models.

Current State 
--------------

Algorithm and visualization codes for **latitude-longitude contour maps**, 
**polar contour maps**, **pressure-latitude zonal mean contour plots**, 
**zonal mean line plots**, and **Cloud Top Height-Tau** joint histograms 
from COSP cloud simulator output. Plots can be created for annual
and seasonal climatologies.

The package also supports custom user diagnostics, by specifying
plot type, desired region (global, ocean, land, etc.), 
pressure levels for variables with the vertical dimension.

For flexibility, the code structure cleanly separates data manipulation 
(reading input files, processing data, etc) from plotting functions.
To satisfy specific user tastes, two graphical back-ends are available: 

* `matplotlib <https://matplotlib.org>`_/ `cartopy <http://scitools.org.uk/cartopy>`_ (**mpl**)
* `UV-CDAT <https://uvcdat.llnl.gov/index.html>`_ VCS (**vcs**)

Additional back-ends could be implemented if the need arose.

+--------------------------------------+---------------------------------------+
| .. figure:: _static/index/fig1.png   | .. figure:: _static/index/fig2.png    |
|    :align: center                    |    :align: center                     |
|    :target: _static/index/fig1.png   |    :target: _static/index/fig2.png    |
|                                      |                                       |
|    Latitude-longitude contour map    |    Polar contour map (mpl)            |
|    (mpl)                             |                                       |
+--------------------------------------+---------------------------------------+
| .. figure:: _static/index/fig3.png   | .. figure:: _static/index/fig4.png    |
|    :align: center                    |    :align: center                     |
|    :target: _static/index/fig3.png   |    :target: _static/index/fig4.png    |
|                                      |                                       |
|    Pressure-latitude zonal mean      |    Zonal mean line plot (vcs)         |
|    contour plot (vcs)                |                                       |
+--------------------------------------+---------------------------------------+
| .. figure:: _static/index/fig5.png                                           |
|    :figwidth: 50 %                                                           |
|    :align: center                                                            |
|    :target: _static/index/fig5.png                                           |
|                                                                              |
|    Cloud Top Height-Tau joint                                                |
|    histograms (mpl)                                                          |
+--------------------------------------+---------------------------------------+

Feature availability for each backend
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. |check| unicode:: U+2714  .. checkmark symbol
.. |ballot| unicode:: U+2718  .. ballot symbol

Not all plot sets and feature are currently supported for every backend.
The table below summarizes current status.

+--------------------------------------------+---------+---------------------+
| Plot set or Feature                        | mpl     | vcs                 |
+============================================+=========+=====================+
| Latitude-longitude contour maps            | |check| | |check|             |
+--------------------------------------------+---------+---------------------+
| Polar contour maps                         | |check| | |check|             |
+--------------------------------------------+---------+---------------------+
| Pressure-latitude zonal mean contour plots | |check| | |check|             |
+--------------------------------------------+---------+---------------------+
| Zonal mean line plots                      | |check| | |check|             |
+--------------------------------------------+---------+---------------------+
| Cloud Top Height-Tau joint histograms      | |check| | |ballot| :sup:`[1]` |
+--------------------------------------------+---------+---------------------+
| Multi-processing                           | |check| | |ballot| :sup:`[2]` |
+--------------------------------------------+---------+---------------------+

| :sup:`[1]` Defaults to mpl instead.
| :sup:`[2]` Your mileage will vary (`Github issue #88 <https://github.com/ACME-Climate/acme_diags/issues/88>`_)

