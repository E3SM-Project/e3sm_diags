
Running ACME Diags on NERSC Edison
==================================

Installation
------------

1. Log on to NERSC Edison

2. Load anaconda module

::

    module load python/2.7-anaconda-4.4

3a. Remove any cached Anaconda downloaded packages. This guarantees you
get the latest packages.

::

    conda clean --all

3b. We'll create an Anaconda environment named ``acme_diags_env`` and
install ``acme_diags``. \* In case you're curious, the command below
installs ``acme_diags`` and all it's dependencies by looking through the
``acme``, ``conda-forge`` (default channel for all software), and
``uvcdat`` channels in this order.

::

    conda create -n acme_diags_env -c acme -c conda-forge -c uvcdat acme_diags

4. Activate the newly created Anaconda environment.

::

    source activate acme_diags_env

5. Install this add-on for vcs.

::

    conda install mesalib -c conda-forge -c uvcdat

-  This is needed because we don't want to use the X windowing system.


Running a simple test
---------------------

Create a parameters file called ``myparams.py``.

::

    touch myparams.py

Copy and paste the below code into ``myparams.py`` using your
favorite text editor. Adjust any options as you like.

.. code:: python

    reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
    test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    sets = ["lat_lon"]

    # optional settings below
    diff_title = 'Model - Obs'

    backend = 'mpl'  # 'mpl' is for the matplotlib plots.

    results_dir = 'lat_lon_demo'  # name of folder where all results will be stored

By default, all of the ACME diagnostics are run for the ``sets`` that
we defined above. This takes some time, so instead we create our own
diagnostics to be run.

Copy and paste the code below in ``mydiags.cfg``. View
`this <./available-parameters.ipynb>`__ document for all available
parameters.

::

    [Diags]
    case_id = "GPCP_v2.2"
    variables = ["PRECT"]
    ref_name = "GPCP_v2.2"
    reference_name = "GPCP (yrs1979-2014)"
    seasons = ["ANN", "DJF"]
    regions = ["global"]
    test_colormap = "WhiteBlueGreenYellowRed.rgb"
    reference_colormap = "WhiteBlueGreenYellowRed.rgb"
    diff_colormap = "BrBG"
    contour_levels = [0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 14, 15, 16]
    diff_levels = [-5, -4, -3, -2, -1, -0.5, 0.5, 1, 2, 3, 4, 5]

    [Diags 2]
    case_id = "SST_CL_HadISST"
    variables = ["SST"]
    ref_name = "HadISST_CL"
    reference_name = "HadISST/OI.v2 (Climatology) 1982-2001"
    seasons = ["ANN", "MAM"]
    contour_levels = [-1, 0, 1, 3, 6, 9, 12, 15, 18, 20, 22, 24, 26, 28, 29]
    diff_levels = [-5, -4, -3, -2, -1, -0.5, -0.2, 0.2, 0.5, 1, 2, 3, 4, 5]

8. Run ACME diags.

::

    acme_diags_driver.py -p myparams.py -d mydiags.cfg

9. Open the following webpage to view the results:

::

    lat_lon_demo/viewer/index.html


Running a full diagnostics suite
--------------------------------

Copy and paste the following into ``myparams.py`` using your
favorite text editor:

.. code:: python

  reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
  test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

  test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

  sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']

  # optional settings below
  diff_title = 'Model - Obs'

  backend = 'mpl'  # 'mpl' is for the matplotlib plots.

  results_dir = 'diag_demo'  # name of folder where all results will be stored

  multiprocessing = True
  num_workers =  24

Compared to the previous short test above, note the following changes:

* Generate plots for all the available sets ('zonal_mean_xy', 'zonal_mean_2d', 
  'lat_lon', 'polar', 'cosp_histogram').
* Turn on multiprocessing with 24 workers.

Since the example above turns on multiprocessing, it should not be run interactively
on the Edison login nodes (NERSC would likely flag it and let you know about it).
Instead, it can be run either in an interactive session on compute nodes, or as a batch
job.


Interactive session on compute nodes
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First, request an interactive session with a single node (24 cores) for one hour
(running this example should take much less than this): ::

  salloc -A acme --nodes=1 --partition=regular --time=01:00:00

Once the session is available, launch ACME Diags: ::

  source activate acme_diags_env
  acme_diags_driver.py -p myparams.py

Batch job
^^^^^^^^^

Alternatively, you can also create a script and submit it to the batch system.
Copy and paste the code below into a file named ``diags.bash``:

.. code:: bash
 
  #!/bin/bash -l
  #SBATCH --job-name=diags
  #SBATCH --output=diags.o%j
  #SBATCH --partition=regular
  #SBATCH --account=acme
  #SBATCH --nodes=1
  #SBATCH --time=01:00:00

  source activate acme_diags_env
  cd /global/cscratch1/sd/golaz/tmp
  acme_diags_driver.py -p myparams.py

And then submit it ::

  sbatch diags.bash

That's it!

