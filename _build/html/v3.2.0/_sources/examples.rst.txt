********
Examples
********

Introduction
============

The model and observation data are located at NERSC, so you can run the examples on Perlmutter CPU.

Make sure you're using version 2.0.0 or greater of e3sm_diags.

Enter a ``conda`` environment that has ``e3sm_diags`` installed.
Follow either `a quickstart guide <quickguides/index.html>`__ or `the instructions here <install.html>`__.

The following tables visualize which comparison each of the following examples demonstrate:

.. list-table:: Time Series Comparison Examples
    :header-rows: 1

    * -
      - Model
      - Observation

    * - Model
      - Ex.1, Ex.2
      - Ex.3

    * - Observation
      -
      -

.. list-table:: Climatology Comparison Examples
    :header-rows: 1

    * -
      - Model
      - Observation

    * - Model
      - Ex.4
      - Ex.5, Ex.6

    * - Observation
      -
      - Ex.7


Examples
========


1. Model Time-series vs Model Time-series: Historical H1 (2011-2013) vs Historical H1 (1850-1852)
-------------------------------------------------------------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex1-model_ts-vs-model_ts>`__ covers how to compare different time slices between two models.
In this case, we're comparing two different three-year time slices from the same model.
The raw model output was run through NCO, which created the time-series files.

2. Model Time-series vs Model Time-series with CMIP data
--------------------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex2-model_ts-vs-model_ts-cmip>`__ covers how to compare different time slices between two models with CMIP5 conventions.
In this case, we're comparing two different three-year time slices from the same model.
The raw model output was run through NCO, which created the time-series files.

3. Model Time-series vs Observation Time-series with CMIP data
--------------------------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex3-model_ts-vs-obs_ts-cmip>`__ covers how to compare different time slices between observational data and a model with CMIP5 conventions.
In this case, we're comparing two different three-year time slices.
The raw model output was run through NCO, which created the time-series files.

4. Model Climatology vs Model Climatology
-----------------------------------------

`This example <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex4-model-vs-model>`__ covers how to compare the climatology between two different model outputs.
The raw model output was run through NCO, which computed the climatology.
We are comparing two simulations: F1850COSP and FC5COSP.

5. Model Climatology vs Observation Climatology
-----------------------------------------------

`This example <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/ex5-model-vs-obs>`__ covers how to compare the climatology between model output data and observational data.
The raw model output was run through NCO, which computed the climatology.
We are comparing model and reanalysis data for surface air temperature for two areas: over land and globally.

6. Model Climatology vs Observation Climatology -- Zonal Mean 2D and Lat/Lon
----------------------------------------------------------------------------

`This example <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex6-model-vs-obs-custom>`__ covers how to compare the climatology between model output data and observational data
on two different sets: ``zonal_mean_2d`` an ``lat_lon``.

7. Observation Climatology vs Observation Climatology
-----------------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex7-obs-vs-obs>`__ covers how to compare observational data with itself,
so you can compare different version of the data, or the same variable from different datasets.
We are comparing CERES EBAF TOA version 2.8 and 4.0.

8. Native Grid Visualization (v3.1.0)
--------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex8-native-grid-visualization>`__ demonstrates how to visualize model data on its native grid
(e.g., cubed-sphere, unstructured grids) without regridding to a regular lat-lon grid.
This feature uses UXarray to preserve native grid features and is particularly useful for high-resolution models with complex grid structures.
The example shows model vs model comparison using snapshot analysis on native grids.

**Key features:**

- Visualize data on native grids without regridding
- Preserve native grid features and characteristics
- Support for cubed-sphere and unstructured grids
- Uses UXarray for grid-aware operations

9. Snapshot Analysis for Core Sets (v3.1.0)
--------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex9-snapshot-analysis>`__ demonstrates time slice analysis on core diagnostic sets.
Instead of computing climatological seasonal means, this analyzes individual time steps from model output using index-based time selection.

This is useful for analyzing specific events, comparing model states at particular time points,
or understanding temporal evolution without time averaging. The example shows how to use the ``time_slices`` parameter
on multiple diagnostic sets (lat_lon, zonal_mean_2d, polar, meridional_mean_2d, zonal_mean_2d_stratosphere).

**Key features:**

- Index-based time selection (e.g., time_slices = ["0", "1", "2"])
- Analyze individual time steps without temporal averaging
- Event-based or process-oriented diagnostics
- Works across multiple core diagnostic sets

**Note:** ``time_slices`` and ``seasons`` parameters are mutually exclusive.

10. Precipitation PDF Diagnostics (v3.2.0)
-------------------------------------------
`This example <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/ex10-precip-pdf>`__ demonstrates how to run precipitation PDF (Probability Density Function) diagnostics
comparing model daily precipitation output against observational datasets (GPCP, IMERG).

The precip_pdf diagnostic set calculates frequency and amount PDFs of daily precipitation and supports:
- Multiple reference datasets for comparison
- Seasonal subsetting (ANN, DJF, MAM, JJA, SON)
- Regional analysis (global, TROPICS, CONUS)
- NetCDF caching for offline use

**Key features:**

- Analyzes daily precipitation PDFs on exponentially-spaced bins
- Compares against GPCP 1DD Daily and GPM IMERG observations
- Generates interactive HTML viewer with detailed plots
- Supports both Python script and command-line execution

Running the Examples
====================

Preprocessing the Data (Ex.2,3 only)
------------------------------------
There are two things you must do to prepare CMIP data to be used by ``e3sm_diags``.

1. Since CMIP model output files follow specific file naming conventions,
they must be renamed to follow the E3SM file naming conventions:
``<variable>_<start_yr>01_<end_yr>12.nc``.

* Ex: renaming ``tas_Amon_CESM1-CAM5_historical_r1i2p1_196001-201112.nc`` to ``tas_196001_201112.nc``.

2. All of the variables should be in the same directory.

If you're running with data not used in this example, you must do the above two steps.

Setting the Parameters
----------------------

The parameters file contains information related to the location
of the data, what years to run the diagnostics on, what plots to create, and more.

The configuration file provides information about the diagnostics you are running.
This is used in Ex.4, 5, 7, 8, 9, 10.

Parameters for each example can be found in
`this directory <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples>`__.

Running the Diagnostics
-----------------------
Enter the directory containing the example you want to run.
Edit the parameter file ``prefix`` value so the results will be placed in your web directory.
Use the code below to run the diagnostics.

    .. code::

        # Allocate a node to run an interactive session on. You can also use a batch job.
        salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm
        # Enter the E3SM Unified environment. For Perlmutter CPU, the command to do this is:
        source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
        # Running Ex.1. For examples 4, 5, 7, 8, 9, 10 append ``-d diags.cfg``.
        python ex1.py --multiprocessing --num_workers=32
        # You may need to change permissions on your web directory to see the example output.
        chmod -R 755 <your web directory>

Note: For Ex.7, you shouldn't run the software
with just ``python ex7.py`` (i.e., without a ``.cfg`` file).
The reason is that ``e3sm_diags`` doesn't support the obs vs obs comparison with all of the
default variables. For each of the plot sets, the user needs to make a ``*_obs_vs_obs.cfg`` file in
`this directory <https://github.com/E3SM-Project/e3sm_diags/tree/master/e3sm_diags/driver/default_diags>`__.

Viewing the Results
-------------------
Results from running all the examples can be found `here <https://portal.nersc.gov/cfs/e3sm/forsyth/examples/>`__.
You can navigate to ``https://portal.nersc.gov/cfs/e3sm/forsyth/examples/<example directory>/viewer/`` to
see the viewer for a specific example.

These were generated with the following script:

    .. code::

       #!/bin/bash
       #SBATCH  --job-name=e3sm_diags_examples
       #SBATCH  --account=e3sm
       #SBATCH  --nodes=1
       #SBATCH  --output=/global/homes/f/forsyth/e3sm_diags_examples.o%j
       #SBATCH  --exclusive
       #SBATCH  --time=1:00:00
       #SBATCH  --qos=regular
       #SBATCH  --constraint=haswell

       cd /global/homes/f/forsyth/e3sm_diags/examples
       # Edit `<your directory>` lines:
       # emacs ex1-model_ts-vs-model_ts/ex1.py
       # emacs ex2-model_ts-vs-model_ts-cmip/ex2.py
       # emacs ex3-model_ts-vs-obs_ts-cmip/ex3.py
       # emacs ex4-model-vs-model/ex4.py
       # emacs ex5-model-vs-obs/ex5.py
       # emacs ex6-model-vs-obs-custom/ex6.py
       # emacs ex7-obs-vs-obs/ex7.py
       # emacs ex8-native-grid-visualization/ex8.py
       # emacs ex9-snapshot-analysis/ex9.py
       # emacs ex10-precip-pdf/ex10.py

       source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
       cd ex1-model_ts-vs-model_ts
       python ex1.py --multiprocessing --num_workers=32
       cd ../ex2-model_ts-vs-model_ts-cmip
       python ex2.py --multiprocessing --num_workers=32
       cd ../ex3-model_ts-vs-obs_ts-cmip
       python ex3.py --multiprocessing --num_workers=32
       cd ../ex4-model-vs-model
       python ex4.py --multiprocessing --num_workers=32 -d diags.cfg
       cd ../ex5-model-vs-obs
       python ex5.py --multiprocessing --num_workers=32 -d diags.cfg
       cd ../ex6-model-vs-obs-custom
       python ex6.py 
       cd ../ex7-obs-vs-obs
       python ex7.py -d diags.cfg
       cd ../ex8-native-grid-visualization
       python ex8.py -d diags.cfg
       cd ../ex9-snapshot-analysis
       python ex9.py -d diags.cfg
       cd ../ex10-precip-pdf
       python ex10.py -d diags.cfg --multiprocessing --num_workers=32
       cd ../

       chmod -R 755 /global/cfs/cdirs/e3sm/www/forsyth/examples
