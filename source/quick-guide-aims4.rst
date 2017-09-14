
A Quick Guide on Running ACME Diags on AIMS4 or ACME1
=====================================================

If you don't enjoy or can't read a lot, just follow this quick guide to
run ``acme_diags`` on ``aims4`` or ``acme1``.

1. Log on to ``aims4``:

::

    ssh -Y aims4.llnl.gov

or ``acme1``:

::

    ssh -Y acme1.llnl.gov

2. If you don't have Anaconda installed, follow `this
guide <https://docs.continuum.io/anaconda/install-linux>`__.

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

4a. Install this addon for vcs.

::

    conda install mesalib -c conda-forge -c uvcdat

-  This is needed because we don't want to use the X windowing system.

5. Create a parameters file called ``myparams.py``.

::

    touch myparams.py

6. Copy and paste the below code into ``myparams.py`` using your
favorite text editor. Adjust any options as you like.

.. code:: python

    reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
    test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    sets = ["lat_lon"]

    # optional settings below
    diff_title = 'Model - Obs'

    backend = 'vcs'  # 'mpl' is for the matplotlib plots.

    results_dir = 'lat_lon_demo'  # name of folder where all results will be stored

7. By default, all of the ACME diagnostics are ran for the ``sets`` that
we defined above. This takes some time, so we'll create our own
diagnostics to be ran. Run the command

::

    touch mydiags.cfg

and paste the code below in ``mydiags.cfg``. View
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

8a. Run the diags.

::

    acme_diags_driver.py -p myparams.py -d mydiags.cfg

8b. You can even run all of the Latitude-Longitude contour plots.

::

    acme_diags_driver.py -p myparams.py

9. Open the following webpage to view the results.

::

    firefox --no-remote lat_lon_demo/viewer/index.html &

-  The ``--no-remote`` option runs this instances of Firefox as a new
   process, thus loading the results webpage faster. ``&`` lets Firefox
   run in the background.

More Options
------------

-  You can modify the ``sets`` parameters in ``myparams.py`` to run
   multiple sets. Possible options are:
   ``'zonal_mean_xy', 'zonal_mean_2d', 'lat_lon, 'polar', 'cosp_histogram'``.
   If the ``sets`` parameter is not defined, all of the aforementioned
   sets are ran. Ex:

   .. code:: python

       sets = ['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon', 'polar', 'cosp_histogram']

-  Diagnostics can be ran in parallel with multi-processing. In
   ``myparams.py``, add ``multiprocessing = True`` and set
   ``num_workers`` to the number of workers you want to use. If
   ``num_workers`` is not defined, it will automatically use 4 processors processes by defualt on a machine. Ex:

   .. code:: python

       # myparams.py
       # In addition to your other parameters, include:
       multiprocessing = True
       num_workers = 4

Below figure shows a scalability test running the package for all lat_lon diagostics on ACME1. Courtesy of Sterling Baldwin. 

.. figure:: _static/_static/quick-guide-aims4/performance_test.png
   :alt: Performance_test

   Figure: Performance test running the package with full set: "lat_lon" diagnostics on ACME1
