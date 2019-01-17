Model Time-series vs Model Time-series: Historical H1 (2011-2013) vs Historical H1 (1850-1852)
----------------------------------------------------------------------------------------------

Introduction and prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide covers how to compare different time slices between two models.
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

Setting the parameters
^^^^^^^^^^^^^^^^^^^^^^

Create a parameters file ``myparams.py`` with the contents below. 
The parameters file contains information related to the location 
of the model time-series data, what years to run the diagnostics 
on, what plots to create, and more parameters.

    .. code:: python

        # Location of the data.
        test_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'
        reference_data_path = '/p/user_pub/work/E3SM/1_0/historical_H1/1deg_atm_60-30km_ocean/atmos/129x256/time-series/mon/ens1/v1/'

        # Set this parameter to True.
        # By default, e3sm_diags expects the test data to be climo data.
        test_timeseries_input = True
        # Years to slice the test data, base this off the years in the filenames.
        test_start_yr = '2011'
        test_end_yr = '2013'

        # Set this parameter to True.
        # By default, e3sm_diags expects the ref data to be climo data.
        ref_timeseries_input = True
        # Years to slice the ref data, base this off the years in the filenames.
        ref_start_yr = '1850'
        ref_end_yr = '1852'

        # When running with time-series data, you don't need to specify the name of the data.
        # But you should, otherwise nothing is displayed when the test/ref name is needed.
        short_test_name = 'historical_H1'
        short_ref_name = 'historical_H1'

        # This parameter modifies the software to accommodate model vs model runs.
        # The default setting for run_type is 'model_vs_obs'.
        run_type = 'model_vs_model'
        # Name of the folder where the results are stored.
        results_dir = 'modTS_vs_modTS_3years'

        # Below are more optional arguments.

        # What plotsets to run the diags on.
        # If not defined, then all available sets are used. 
        sets = ['lat_lon']
        # What seasons to run the diags on.
        # If not defined, diags are ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
        seasons = ['ANN']
        # Title of the difference plots.
        diff_title = 'Model (2011-2013) - Model (1850-1852)'
        # For running with multiprocessing.
        multiprocessing = True
        num_workers = 32

Running the diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Use the command below to run the diagnostics.
Again, if you want to run the container, read the 'Introduction and prerequisites'
section above for instructions on how to do so.

    .. code::

        e3sm_diags -p myparams.py


Viewing the results
^^^^^^^^^^^^^^^^^^^

The results of this example are `on the e3sm_diags documentation site
<../../../sample_output/modTS_vs_modTS_3years/viewer/index.html>`__.

