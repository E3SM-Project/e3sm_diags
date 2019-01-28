Model Time-series vs Model Time-series with CMIP data
-----------------------------------------------------

Introduction and prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide covers how to compare different time slices between two models with CMIP5 conventions.
In this case, we're comparing two different three-year time slices from the same model.
The raw model output was ran through NCO, which created the time-series files.
The model data is located at NERSC, so you can run this example on Cori or Edison.

First, make sure you're using version 1.6.0 or greater of e3sm_diags.

Then make sure you're either:

* In an environment with e3sm_diags installed.
   * Either follow `a quickstart guide <../quickguides/index.html>`__
     or `the instructions here <../install.html>`__
* Or have the container downloaded, and download `this script <https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py>`__ to run it.
   * If you're using containers, when you run the software, replace every instance of
     'e3sm_diags' in this guide with ``python e3sm_diags_container.py --<your_container_runtime>``.
   * See `this guide <../quickguides/quick-guide-cori.html>`__ for more information.


Preprocessing the data
^^^^^^^^^^^^^^^^^^^^^^
There are two things you must do to prepare CMIP data to be used by e3sm_diags.

1. Since CMIP model output files follow specific file naming conventions,
they must be renamed to follow the E3SM file naming conventions, where you have
``<variable>_<start_yr>01_<end_yr>12.nc``.

* Ex: renaming ``tas_Amon_CESM1-CAM5_historical_r1i2p1_196001-201112.nc`` to ``tas_196001_201112.nc``.

2. All of the variables should be in the same directory.

If you're running with data not used in this example,
you must do the above two steps.


Setting the parameters
^^^^^^^^^^^^^^^^^^^^^^

Create a parameters file ``myparams.py`` with the contents below. 
The parameters file contains information related to the location 
of the model time-series data, what years to run the diagnostics 
on, what plots to create, and more parameters.

    .. code:: python

        # Location of the data.
        test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/time-series/CESM1-CAM5_cmip/'
        reference_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/time-series/CESM1-CAM5_cmip'

        # Set this parameter to True.
        # By default, e3sm_diags expects the test data to be climo data.
        test_timeseries_input = True
        # Years to slice the test data, base this off the years in the filenames.
        test_start_yr = '2001'
        test_end_yr = '2003'

        # Set this parameter to True.
        # By default, e3sm_diags expects the ref data to be climo data.
        ref_timeseries_input = True
        # Years to slice the ref data, base this off the years in the filenames.
        ref_start_yr = '1850'
        ref_end_yr = '1852'

        # When running with time-series data, you don't need to specify the name of the data.
        # But you should, otherwise nothing is displayed when the test/ref name is needed.
        short_test_name = 'CESM1-CAM5-historical'
        short_ref_name = 'CESM1-CAM5-historical'

        # This parameter modifies the software to accommodate model vs model runs.
        # The default setting for run_type is 'model_vs_obs'.
        run_type = 'model_vs_model'
        # Name of the folder where the results are stored.
        results_dir = 'modTS_vs_modTS_CMIP_3years'

        # Below are more optional arguments.

        # What plotsets to run the diags on.
        # If not defined, then all available sets are used. 
        sets = ['lat_lon']
        # What seasons to run the diags on.
        # If not defined, diags are ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
        seasons = ['ANN']
        # Title of the difference plots.
        diff_title = 'Model (2001-2003) - Model (1850-1852)'
        ## For running with multiprocessing.
        #multiprocessing = True
        #num_workers = 32

Running the diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Use the command below to run the diagnostics.
Again, if you want to run the container, read the 'Introduction and prerequisites'
section above for instructions on how to do so.

    .. code::

        e3sm_diags -p myparams.py

This run includes all variables predefined in the default configuration file. To make it run fast, please use multiprocessing and run it either in an interactive session on compute nodes, or as a batch job.

