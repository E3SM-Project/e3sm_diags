Observation Climatology vs Observation Climatology
--------------------------------------------------

Introduction and prerequisites
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This guide covers how to compare observational data with itself.
So you can compare different version of the data, or the same variable from different datasets.
We are comparing CERES EBAF TOA version 2.8 and 4.0.
The observational data is located at NERSC, so you can run this example on Cori or Edison.

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
of the obs data, what plots to create, and more parameters.

    .. code:: python
    
        # Location of the ref data.
        reference_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/'
        # Name of the ref obs data, used to find the climo files.
        ref_name = 'ceres_ebaf_toa_v2.8'

        # Location of the test data.
        test_data_path = '/global/cfs/cdirs/e3sm/acme_diags/obs_for_e3sm_diags/climatology/'
        # Name of the test obs data, used to find the climo files.
        test_name = 'ceres_ebaf_toa_v4.0'

        # Name of the folder where the results are stored.
        results_dir = 'obs_vs_obs'
        # What plotsets to run the diags on.
        sets = ['lat_lon']

        # Below are more optional arguments.

        # What seasons to run the diags on.
        # If not defined, diags is ran on ['ANN', 'DJF', 'MAM', 'JJA', 'SON'].
        seasons = ['ANN']
        # 'mpl' is to create matplotlib plots, 'vcs' is for vcs plots.
        backend = 'mpl'


The ``mydiags.cfg`` below provides information about the diagnostics you are running.

    .. code:: ini

        [#]
        sets = ["lat_lon"]
        case_id = "lat_lon_obs_vs_obs"
        ref_name = "ceres_ebaf_toa_v2.8"
        reference_name = "CERES-EBAF"
        variables = ["SWCF"]
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [-180, -160, -140, -120, -100, -80, -60, -40, -20,  0]
        diff_levels = [-60, -50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50, 60]

        [#]
        sets = ["lat_lon"]
        case_id = "lat_lon_obs_vs_obs"
        ref_name = "ceres_ebaf_toa_v2.8"
        reference_name = "CERES-EBAF"
        variables = ["LWCF"]
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [0, 10, 20, 30, 40, 50, 60, 70, 80]
        diff_levels = [-35, -30, -25, -20, -15, -10, -5, -2, 2, 5, 10, 15, 20, 25, 30, 35]

        [#]
        sets = ["lat_lon"]
        case_id = "lat_lon_obs_vs_obs"
        ref_name = "ceres_ebaf_toa_v2.8"
        reference_name = "CERES-EBAF"
        variables = ["NETCF"]
        seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]
        contour_levels = [-135, -120, -105, -90, -75, -60, -45, -30, -15, 0, 15, 30, 45]
        diff_levels = [-75, -50, -40, -30, -20, -10, -5, 5, 10, 20, 30, 40, 50, 75]



Running the diagnostics
^^^^^^^^^^^^^^^^^^^^^^^

Use the command below to run the package with the diags you've created in ``mydiags.cfg``.
Again, if you want to run the container, read the 'Introduction and prerequisites'
section above for instructions on how to do so.

    .. code::

        e3sm_diags -p myparams.py -d mydiags.cfg


Note: Unlike the other examples, you shouldn't run the software
with ``e3sm_diags -p myparams.py``. The reason for this is that
e3sm_diags doesn't support the obs vs obs comparison with all of the
default variables. For each of the plotsets, user needs to make a
``*_obs_vs_obs.cfg`` file in
`this directory <https://github.com/E3SM-Project/e3sm_diags/tree/master/acme_diags/driver/default_diags>`__.
