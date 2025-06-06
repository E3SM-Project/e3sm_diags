..
    Comment: If you want to edit `quick-guide-{machine_name}.rst`, edit `quick-guide-generic.rst` instead and run `generate_quick_guides.py`.

Perlmutter quick guide for running e3sm_diags v3
=========================================================================


1. Installation
-----------------------------------------------------------

We will use the e3sm_unifed environment to install.
For the latest stable release or if you don't have access to e3sm analysis machines,
please instead refer to :ref:`Latest stable release <install_latest>`.

Most of the E3SM analysis software is maintained with an Anaconda metapackage
(E3SM unified environment).
If you have an account on Perlmutter,
then to get all of the tools in the metapackage in your path,
use the activation command below.
(Change ``.sh`` to ``.csh`` for csh shells.)

Below, we also provide the paths for observational data needed by ``e3sm_diags`` (<obs_path>),
and some sample model data for testing (<test_data_path>).
Both <obs_path> and <test_data_path> have two subdirectories:
``/climatology`` and ``/time-series`` for climatology and time-series data respectively.

Also listed below are paths where the HTML files (<html_path>) must be located to be displayed
at their corresponding web addresses (<web_address>).

<activation_command>: ``source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh``

<obs_path>: ``/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/``

<test_data_path>: ``/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/``

<html_path>: ``/global/cfs/cdirs/e3sm/www/<username>/``

<web_address>: ``http://portal.nersc.gov/cfs/e3sm/<username>/``
     


2. Config and run
--------------------------------------------------------

.. _Perlmutter_lat_lon:

Running the annual mean latitude-longitude contour set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Copy and paste the below code into ``run_e3sm_diags.py`` using your favorite text editor.
Adjust any options as you like.

   **Tip:** Some of E3SM's analysis machines (**Acme1, Anvil, Compy, Cori**)
   have web servers setup to host html results.
   On Perlmutter,
   create the directory ``/global/cfs/cdirs/e3sm/www/<username>/`` using your username.
   Set ``results_dir`` to ``/global/cfs/cdirs/e3sm/www/<username>/doc_examples/lat_lon_demo``
   in ``run_e3sm_diags.py`` below. Then, you can view results via a web browser here:
   http://portal.nersc.gov/cfs/e3sm/<username>/doc_examples/lat_lon_demo


    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'
        param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]   #all seasons ["ANN","DJF", "MAM", "JJA", "SON"] will run,if comment out"

        prefix = '/global/cfs/cdirs/e3sm/www/<username>/doc_examples/'
        param.results_dir = os.path.join(prefix, 'lat_lon_demo')
        # Use the following if running in parallel:
        #param.multiprocessing = True
        #param.num_workers = 32
        
        # Use below to run all core sets of diags:
        #runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']
        # Use below to run lat_lon map only:
        runner.sets_to_run = ['lat_lon']
        runner.run_diags([param])


Run in serial with:

    ::

        python run_e3sm_diags.py

The above run has the same results as running ``e3sm_diags -p lat_lon_params.py``
using the code below for ``lat_lon_params.py``:


    .. code:: python

        reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'
        test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/climatology/'

        test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

        sets = ["lat_lon"]
        seasons = ["ANN"]

        # Name of folder where all results will be stored.
        results_dir = '/global/cfs/cdirs/e3sm/www/<username>/doc_examples/lat_lon_demo'

The new way of running (no ``-p``) is implemented in version 2.0.0,
preparing ``e3sm_diags`` to accomodate more diagnostics sets with set-specific parameters.


To enable multiprocessing rather than running in serial, the program will need to be run in an
**interactive session** on compute nodes, or as a **batch job**. 


Interactive session on compute nodes
'''''''''''''''''''''''''''''''''''''

First, request an interactive session with a single node
(128 cores with Perlmutter CPU)
for one hour (running this example should take much less than this).
If obtaining a session takes too long, try to use the ``debug`` partition.
Note that the maximum time allowed for that partition is ``00:30:00``.

    ::

        salloc --nodes 1 --qos interactive --time 01:00:00 --constraint cpu --account=e3sm


Once the session is available, launch E3SM Diagnostics, to activate ``e3sm_unified``:

    ::

        source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh
        python run_e3sm_diags.py --multiprocessing --num_workers=128


We could have also set these multiprocessing parameters in the ``run_e3sm_diags.py`` as well
but we're showing that you can still submit parameters via the command line.

Batch job
'''''''''

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``.

    .. code:: bash
    
        #!/bin/bash -l
        #SBATCH --job-name=diags
        #SBATCH --output=diags.o%j
        #SBATCH --partition=regular
        #SBATCH --account=e3sm
        #SBATCH --nodes=1
        #SBATCH --time=01:00:00
        #SBATCH -C cpu

        source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh 
        python run_e3sm_diags.py --multiprocessing --num_workers=32

And then submit it:

    ::

        sbatch diags.bash

View the status of your job with ``squeue -u <username>``.
Here's the meaning of some values under the State (``ST``) column:

* ``PD``: Pending
* ``R``: Running
* ``CA``: Cancelled
* ``CD``: Completed
* ``F``: Failed
* ``TO``: Timeout
* ``NF``: Node Failure

View results on the web
'''''''''''''''''''''''
Once the run is completed,
open  ``http://portal.nersc.gov/cfs/e3sm/<username>/doc_examples/lat_lon_demo/viewer/index.html`` to view the results.
If you don't see the results, you may need to set proper permissions.
Run ``chmod -R 755 /global/cfs/cdirs/e3sm/www/<username>/``.

**Tip:** Once you're on the webpage for a specific plot, click on the
'Output Metadata' drop down menu to view the metadata for the displayed plot.
Running that command allows the displayed plot to be recreated.
Changing any of the options will modify just that resulting figure.



Running all the core diagnostics sets
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Core diagnostics set includes:
**lat_lon**, **zonal_mean_xy**, **zonal_mean_2d**, **polar**, **cosp_histogram**,
**meridional_mean_2d**.
These diags share a common parameter space (core parameters).
To run all these sets without defining set-specific parameters
(e.g. **plev** for **zonal_mean_2d** and **meridional_mean_2d**.),
replace the ``runner.sets_to_run`` line in ``run_e3sm_diags.py`` with the one below:

 ::

   runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d']


Running area mean time series set
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In v2.0.0, the time series set was implemented to support regional averaged time series plotting
using monthly mean time series input.
This set is enabled if monthly mean time series is processed as documented
:doc:`here <../input-data-requirement>`.

A ``run_e3sm_diags.py`` example for running area mean time series alone:

    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from e3sm_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/'
        param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/time-series/E3SM_v1/'
        param.test_name = 'e3sm_v1'
        
        prefix = '/global/cfs/cdirs/e3sm/www/<username>/doc_examples/'
        param.results_dir = os.path.join(prefix, 'area_mean_with_obs')
        # Use the following if running in parallel:
        #param.multiprocessing = True
        #param.num_workers =  40
        
        # We're passing in this new object as well, in
        # addition to the CoreParameter object.
        
        ts_param = AreaMeanTimeSeriesParameter()
        #ts_param.ref_names = ['none']   # Using this setting will plot only the model data, not the observation data
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['area_mean_time_series']
        runner.run_diags([param, ts_param])


This set can also be ran with the core diagnostics sets,
so that all the plots are shown in one viewer.
The following is an example to run all sets:

    .. code:: python

        import os
        from e3sm_diags.parameter.core_parameter import CoreParameter
        from e3sm_diags.parameter.area_mean_time_series_parameter import AreaMeanTimeSeriesParameter
        from e3sm_diags.run import runner
        
        param = CoreParameter()
        
        param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/climatology/'
        param.test_data_path = '/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.multiprocessing = True
        param.num_workers = 40
        prefix = '/global/cfs/cdirs/e3sm/www/<username>/doc_examples'
        param.results_dir = os.path.join(prefix, 'all_sets')
        
        #
        ##Set specific parameters for new sets
        ts_param = AreaMeanTimeSeriesParameter()
        ts_param.reference_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/'
        ts_param.test_data_path = '/global/cfs/cdirs/e3sm/diagnostics/observations/Atm/time-series/E3SM_v1/'
        ts_param.test_name = 'e3sm_v1'
        ts_param.start_yr = '2002'
        ts_param.end_yr = '2008'
        
        runner.sets_to_run = ['lat_lon','zonal_mean_xy', 'zonal_mean_2d', 'polar', 'cosp_histogram', 'meridional_mean_2d', 'area_mean_time_series']
        runner.run_diags([param, ts_param])


Advanced: Running custom diagnostics
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The following steps are for 'advanced' users, who want to run custom diagnostics.
So, most users will not run the software like this.


By default, with ``e3sm_diags``,
a built in set of variables are defined for each diagonostics sets.
To do a short run, e.g. only running through a subset of variables,
a configuration file is needed to customize the run.


In the following example,
only precipitation and surface sea temperature are run to compare with
model and obs for lat_lon set.
Create ``mydiags.cfg`` file as below.

Check :doc:`Available Parameters <../available-parameters>` for all available parameters.

For a larger configuration file example, look
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/e3sm_diags/driver/default_diags/lat_lon_model_vs_obs.cfg>`_
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


Run E3SM diagnostics with the ``-d`` parameter.
Use the :ref:`above run script <Perlmutter_lat_lon>`. And run as following:

    ::

        python run_e3sm_diags.py -d mydiags.cfg


