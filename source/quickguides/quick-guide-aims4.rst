
Acme1 quick guide for running e3sm_diags v2
=========================================================================

Log on to ``acme1``:

::

    ssh -Y acme1.llnl.gov



1. Installation
-----------------------------------------------------------

For now, we recommend two methods to install:



1a. Installation via conda
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you want to install the latest stable version of ``e3sm_diags``, please follow :ref:`Latest stable release <install_latest>` to install via ``conda``. Remember to install conda/miniconda or load the anaconda module of the machine

Make sure you're using ``bash``
 


1b. Installation: via e3sm_unified environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Most of the E3SM analysis software is maintained with an Anaconda metapackage. If you have an account on an E3SM supported machines(**Cori, Compy, Acme1, Anvil, Cooley, Rhea**), to get all of the tools in the metapackage in your path, use the command below and the activation path from below for Acme1. 


**acme1**
    ::

     source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified.sh

<obs_path> at:``/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags``

<test_data_path> at:``/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags``

<html_path> at: ``/var/www/acme/acme-diags/<username>``

<web_address> at: ``https://acme-viewer.llnl.gov/<username>``



2. Config and run
--------------------------------------------------------

Running the annual mean latitude-longitude contour set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy and paste the below code into ``run_e3sm_diags.py`` using your favorite text editor. Adjust any options as you like.

   **Tip:** Some of E3SM's analysis machines (**Cori, Compy, Acme1, Anvil**) have web servers setup to host html results. For instance, on Acme1, make a folder in the following directory ``/var/www/acme/acme-diags/`` based off your username.
   Then you can set ``results_dir`` to  ``/var/www/acme/acme-diags/<username>/lat_lon_demo`` in ``run_e3sm_diags.py`` below
   to view the results via a web browser here: http://portal.nersc.gov/project/acme/<username>/lat_lon_demo


    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]   #all seasons ["ANN","DJF", "MAM", "JJA", "SON"] will run,if comment out"

        prefix = '/var/www/acme/acme-diags/zhang40/tests/'
        param.results_dir = os.path.join(prefix, 'lat_lon_demo')
        param.multiprocessing = True
        param.num_workers = 32
        
        #use below to run all core sets of diags:
        #runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
        #use below to run lat_lon map only
        runner.sets_to_run = ['lat_lon']
        runner.run_diags([param])


Run by following:

    ::

        python run_e3sm_diags.py

This new way of running is implemented in version 2.0.0, in order to prepare ``e3sm_diags`` for accomodating more diagnostics sets with set-specific parameters. The above run has the same results has :ref:`the parameter file linked here <cori-params-v1>`, which was run using ``e3sm_diags -p lat_lon_demo.py``.

Once the run is completed, open the following webpage served at : https://acme-viewer.llnl.gov/<username> to view the results.


    ::

        lat_lon_demo/viewer/index.html

**Tip:** Once you're on the webpage for a specific plot, click on the
'Output Metadata' drop down menu to view the metadata for the displayed plot.
Running that command allows the displayed plot to be recreated.
Changing any of the options will modify the just that resulting figure.



Running all the core diagnostics sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Core diagnostics set includes: **lat_lon**, **zonal_mean_xy**, **zonal_mean_2d**, **polar**, **cosp_histogram**, **meridional_mean_2d**. These diags share a common parameter space (core parameters). To run all these sets without defining set-specific parameters (i.e. **plev** for **zonal_mean_2d** and **meridional_mean_2d**.), use following in ``run_e3sm_diags.py``:

 ::

   runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']


Running area mean time series set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In v2.0.0, the time series set was implemented to support regional averaged time series plot using monthly mean time series input. This set is enabled if monthly mean time series is processed as documented :doc:`here <../input-data-requirement>`.

A ``run_e3sm_diags.py`` example for running area mean time series alone:

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from acme_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/time-series/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/time-series/E3SM_v1/'
        param.test_name = 'e3sm_v1'
        
        prefix = '/var/www/acme/acme-diags/zhang40/tests/'
        param.results_dir = os.path.join(prefix, 'area_mean_with_obs')
        #param.multiprocessing = True
        #param.num_workers =  40
        
        # We're passing in this new object as well, in
        # addition to the CoreParameter object.
        
        ts_param = AreaMeanTimeSeriesParameter()
        #ts_param.ref_names = ['none']   #This setting plot model data only
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['area_mean_time_series']
        runner.run_diags([param, ts_param])


This set can also be ran with the core diagnostics sets, so that all the plots are shown in one viewer. Following is an example to run all sets:

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from acme_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.multiprocessing = True
        param.num_workers = 40
        prefix = '/global/project/projectdirs/acme/www/zhang40'
        param.results_dir = os.path.join(prefix, 'all_sets')
        
        #
        ##Set specific parameters for new sets
        ts_param = AreaMeanTimeSeriesParameter()
        ts_param.reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/time-series/'
        ts_param.test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/time-series/E3SM_v1/'
        ts_param.test_name = 'e3sm_v1'
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'area_mean_time_series']


Advanced: Running custom diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following steps are for 'advanced' users, who want to run custom diagnostics.
So most users will not run the software like this.


By default, with ``e3sm_diags``, a built in set of variables are defined for each diagonostics sets. To do a short run, i.e. only run through a subset of variables, the a configuration files is needed to customize the run.


In the following example, only precipitation and surface sea temperature are ran to compare with model and obs for lat_lon set. Create a ``mydiags.cfg`` file as following.

Check :doc:`Available Parameters <../available-parameters>`
for all available parameters.

For more examples of these types of files, look
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/default_diags/lat_lon_model_vs_obs.cfg>`_
for the cfg file that was used to create all of the latitude-longitude sets.


    ::

        [#]
        sets = ["lat_lon"]
        case_id = "GPCP_v2.3"
        variables = ["PRECT"]
        ref_name = "GPCP_v2.3"
        reference_name = "GPCP"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        regions = ["global"]
        test_colormap = "WhiteBlueGreenYellowRed.rgb"
        reference_colormap = "WhiteBlueGreenYellowRed.rgb"
        diff_colormap = "BrBG"
        contour_levels = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]
        
        
        [#]
        sets = ["lat_lon"]
        case_id = "SST_HadISST"
        variables = ["SST"]
        ref_name = "HadISST"
        reference_name = "HadISST/OI.v2"
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [-1, 0, 1, 3, 6, 9, 12, 15, 18, 20, 22, 24, 26, 28, 29]
        diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]


Run E3SM diagnostics with the ``-d`` parameter. And use the :ref:`run script linked here <cori-params-v2>`. And run as following:

    ::

        python run_e3sm_diags.py -d mydiags.cfg


