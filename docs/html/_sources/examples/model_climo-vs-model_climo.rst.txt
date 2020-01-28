Model Climatology vs Model Climatology
--------------------------------------

Introduction and prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide covers how to compare the climatology between two different model output.
The raw model output was ran through NCO, which computed the climatology.
You can see in the parameters below we are comparing two simulations: F1850COSP and FC5COSP.
The model data is located at NERSC, so you can run this example on Cori or Edison.

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
of the model data, what plots to create, and more parameters.

    .. code:: python

        # Location of the ref data.
        reference_data_path = '/global/cfs/cdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/'
        # Name of the ref model data, used to find the climo files.
        ref_name = '20161118.beta0.F1850COSP.ne30_ne30.edison'
        # An optional, shorter name to be used instead of the ref_name.
        short_ref_name = 'Ref: beta0.F1850COSP_ne30'

        # Location of the test data.
        test_data_path = '/global/cfs/cdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology'
        # Name of the test model data, used to find the climo files.
        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        # An optional, shorter name to be used instead of the test_name.
        short_test_name = 'Test: beta0_FC5COSP_ne30'

        # What plotsets to run the diags on.
        sets = ['lat_lon']
        # Name of the folder where the results are stored.
        results_dir = 'model_to_model'
        # This parameter modifies the software to accommodate model vs model runs.
        # The default setting for run_type is 'model_vs_obs'.
        run_type = 'model_vs_model' 

        # Below are more optional arguments.

        # 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
        backend = 'mpl'
        # Title of the difference plots.
        diff_title = 'Test Model - Ref Model'
        # For running with multiprocessing.
        #multiprocessing = True
        #num_workers = 32

The ``mydiags.cfg`` below provides information about the diagnostics you are running.
We have two runs with two variables (PRECT and SST) with all seasons selected.

    .. code:: ini

        [#]
        sets = ["lat_lon"]
        case_id = "model_vs_model"
        variables = ["PRECT"]
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        regions = ["global"]
        test_colormap = "WhiteBlueGreenYellowRed.rgb"
        reference_colormap = "WhiteBlueGreenYellowRed.rgb"
        diff_colormap = "BrBG"
        contour_levels = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]

        [#]
        sets = ["lat_lon"]
        case_id = "model_vs_model"
        variables = ["SST"]
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [-1, 0, 1, 3, 6, 9, 12, 15, 18, 20, 22, 24, 26, 28, 29]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]



Running the diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Use the command below to run the package with the diags you've created in ``mydiags.cfg``.
Again, if you want to run the container, read the 'Introduction and prerequisites'
section above for instructions on how to do so.

    .. code::

        e3sm_diags -p myparams.py -d mydiags.cfg


To run the package with the complete variable list, please use multiprocessing and run it either in an interactive session on compute nodes, or as a batch job.
