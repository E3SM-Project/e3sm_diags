Configuration and Running 
=========================

Configuration
-------------

Before runnning the diagnostics, a python script must be prepared to specify parameters for running the diags. 


1. Create a Python script, ex: ``myparams.py``. This script contains simply
   pairs of keys and values. **At minimum, you must define values for the following:**

-  **reference_data_path: path to the reference (observational)
   data**
-  **test_data_path: path to the test (model output) data**
-  **test_name: name of the test (model output) file**

2. There are many other parameters that allow the user to customize
   regridding method, plotting backend, and much more. See
   :doc:`available parameters <available-parameters>` for a description.

An example ``myparams.py`` script is shown below.

.. code:: python

    reference_data_path = '/p/user_pub/e3sm/e3sm_diags_data/obs_for_e3sm_diags/climatology/'
    test_data_path = '/p/user_pub/e3sm/e3sm_diags_data/test_model_data_for_acme_diags/climatology/'
    
    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
    
    # a few optional parameters
    #backend = 'vcs'  # use 'mpl' for matplotlib as default
    
    results_dir = 'myresults'  # name of folder where all results will be stored


Beginning version **v2**, we **added** a new way to run the diags using a python script via an API manner. The same pair of keys and values are required but the Python script is slightly different. An example python interface script as follows `run_e3sm_diags.py`:

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


To add your own diagnostics, create a cfg file like the one below, which
we call ``mydiags.cfg``. (All :doc:`available parameters <available-parameters>` 
can also serve as keys in the cfg file.)

.. code::

    [#]
    # What sets to run this diagnostics on
    sets = ['lat_lon']
    
    # Diagnostics results are saved in a folder named after the case_id
    case_id = "lat_lon_MERRA"
    
    # variables, ref_name, and season are keywords for obs file searching 
    variables = ["T"]  
    ref_name = "MERRA"
    seasons = ["ANN", "JJA"]
    
    # Name of the observation that will appear on the output plot
    reference_name = "MERRA Analysis 1979-2013 NASA"
    
    # User-specified pressure levels
    plevs = [200.0, 850.0]
    
    # User-defined regions, the default region is "global" if region is empty
    # Find default_regions.py in this repo for a list of all possible regions
    regions = ["land", "ocean_TROPICS"] 

If you have multiple diagnostics you want to run, create a cfg file with multiple
entries:

.. code::

    [#]
    # put all of the parameters for a diags run here

    [#]
    # another diags run

``[#]`` is a special title can be used for multiple runs. When used as a title, you don't need to create a new, unique
title for that diagnostics run.


Running
-------

If you **don't** have your own diagnostic file (e.g. ``mydiags.cfg``), simply run: ::

  e3sm_diags -p myparams.py


to generate the standard set of E3SM diagnostics figures.
If you do have your own own diagnostic file, specify it on the command line: ::

  e3sm_diags -p myparams.py -d mydiags.cfg

If you have `e3sm_diags` v2 installed, run: ::

  python run_e3sm_diags.py

Or ::

  python run_e3sm_diags.py -d mydiags.cfg

View the results by opening ``index.html`` in the location specified.


Model vs Model Run
~~~~~~~~~~~~~~~~~~

View
`this <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/model-vs-model/model-vs-model.ipynb>`__ Jupyter Note book and the parameter files ``myparams.py`` and ``mydiags.cfg`` residing in the same directory.

Observation vs Observation Run 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
View
`this <https://github.com/E3SM-Project/e3sm_diags/tree/master/examples/obs-vs-obs/obs-vs-obs.ipynb>`__ Jupyter Note book and the parameter files ``myparams.py`` and ``mydiags.cfg`` residing in the same directory.


Model vs Observation Run 
~~~~~~~~~~~~~~~~~~~~~~~~

Providing model versus observation datasets is the basic use of this software, examples can be found in the quick start guides. 
Check `this <https://github.com/E3SM-Project/e3sm_diags/blob/master/examples/model-vs-obs/model-vs-obs.ipynb>`__ Jupyter Note book and the 
parameter files ``myparams.py`` and ``mydiags.cfg`` residing in the same directory. 
