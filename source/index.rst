.. ACME Diags documentation master file, created by
   sphinx-quickstart on Tue Sep  5 17:27:47 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

************************
ACME Diagnostics Package
************************

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   self
   quick-guide-aims4
   quick-guide-edison
   install-config-run
   available-parameters
   add-new-diagnostics
   contributing


Overview
--------

This diagnostics package is constructed to support the diagnostics
needs of DOE's `Accelerated Climate Modeling for Energy (ACME)
project <https://climatemodeling.science.energy.gov/projects/accelerated-climate-modeling-energy>`__.
The ultimate goal of this work is to develop a comprehensive diagnostics package
that:

-  fully integrates the functionality of NCAR's AMWG diagnostics
   package;
-  utilizes most updated observational datasets, including remote
   sensing, reanalysis and in-situ datasets;
-  interfaces with diagnostics developed from different ACME focus
   groups: atmosphere group, coupled simulation group, land group;
-  interacts effectively with the PCMDI's metrics package and the ARM
   diagnostics package through a unifying framework: `Community
   Diagnostics Package (CDP) <https://github.com/UV-CDAT/cdp>`_.
-  is flexible for user-specified diagnostics and configuration for
   use by other climate models.

Current State 
--------------

Algorithm and visualization codes for **latitude-longitude contour maps**, 
polar contour maps, **pressure-latitude zonal mean contour plots**, 
**zonal mean line plots**, and **Cloud Top Height-Tau** joint histograms 
from COSP cloud simulator output. Plots can be created for annual
and seasonal climatologies.

The package also supports custom user diagnostics, by specifying
plot type, desired region (global, ocean, land, etc.), 
pressure levels for variables with the vertical dimension.

For flexibility, the code structure cleanly separates data manipulation 
(reading input files, processing data, etc) from plotting functions.
To satisfy specific user tastes, two graphical back-ends are available: 

* `matplotlib <https://matplotlib.org>`_/ `cartopy <http://scitools.org.uk/cartopy>`_
* `UV-CDAT <https://uvcdat.llnl.gov/index.html>`_ VCS

Additional back-ends could be implemented if the need arose.

+--------------------------------------+---------------------------------------+
| .. figure:: _static/index/fig1.png   | .. figure:: _static/index/fig2.png    |
|    :align: center                    |    :align: center                     |
|    :target: _static/index/fig1.png   |    :target: _static/index/fig2.png    |
|                                      |                                       |
|    Latitude-longitude contour map    |    Polar contour map (vcs)            |
|    (vcs)                             |                                       |
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
|    histograms (matplotlib)                                                   |
+--------------------------------------+---------------------------------------+


