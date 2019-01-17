Model Climatology vs Observation Climatology
--------------------------------------------

Introduction and prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide covers how to compare the climatology between model output data and observational data.
The raw model output was ran through NCO, which computed the climatology.
We are comparing the surface air temperature averaged over land and global.
The model and observational data are located at NERSC, so you can run this example on Cori or Edison.

Make sure you're either:

* In an environment with e3sm_diags installed.
   * Either follow `a quickstart guide <https://e3sm-project.github.io/e3sm_diags/docs/html/quickguides/index.html>`__
     or `the instructions here <https://e3sm-project.github.io/e3sm_diags/docs/html/install.html>`__
* Or have the container downloaded, and download `this script <https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/acme_diags/container/e3sm_diags_container.py>`__ to run it.
   * If you're using containers, when you run the software, replace every instance of
     'e3sm_diags' in this guide with ``python e3sm_diags_container.py --<your_container_runtime>``.
   * See `this guide <../quickguides/quick-guide-cori.html>`__ for more information.

Setting the parameters
^^^^^^^^^^^^^^^^^^^^^^

Create a parameters file ``myparams.py`` with the contents below. 
The parameters file contains information related to the location 
of the model and obs data, what plots to create, and more parameters.

    .. code:: python

        # Location of the data.
        reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
        test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
        # Name of the test model data, used to find the climo files.
        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        # An optional, shorter name to be used instead of the test_name.
        short_test_name = 'beta0.FC5COSP.ne30'

        # What plotsets to run the diags on.
        sets = ['lat_lon']
        # Name of the folder where the results are stored.
        results_dir = 'era_tas_land'

        # Below are more optional arguments.

        # 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
        backend = 'mpl'
        # Title of the difference plots.
        diff_title = 'Model - Obs.'
        # Save the netcdf files for each of the ref, test, and diff plot.
        save_netcdf = True
        # For running with multiprocessing.
        multiprocessing = True
        num_workers = 32


The ``mydiags.cfg`` below provides information about the diagnostics you are running.

    .. code:: ini

        [#]
        sets = ["lat_lon"]
        case_id = "ERA-Interim"
        variables = ["TREFHT"]
        regions = ["land", "global"]
        ref_name = "ERA-Interim"
        reference_name = "ERA-Interim Reanalysis 1980-2016"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [-35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25, 30, 35, 40]
        diff_levels = [-15, -10, -5, -2, -1, -0.5, -0.2, 0, 0.2, 0.5, 1, 2, 5, 10, 15]



Running the diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Use the command below to run the package with the diags you've created in ``mydiags.cfg``.
Again, if you want to run the container, read the 'Introduction and prerequisites'
section above for instructions on how to do so.

    .. code::

        e3sm_diags -p myparams.py -d mydiags.cfg


To run the package with the complete variable list, use the command below.

    .. code::

        e3sm_diags -p myparams.py
