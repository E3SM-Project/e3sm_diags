Running via the API
===================

As of v2.0.0 of ``e3sm_diags``, the new way to run the software is via the API.
The reason for this is because of the new set-specific parameters being added
to the new plotsets.

However, there is backwards compatibility, so you can still run ``e3sm_diags``
as you normally would in v2.0.0.
Only the current plot sets will be supported running the "old way", since
the newer plotsets, like ``area_mean_time_series``, have set-specific parameters
and you need to run it via the API.

Getting the Environment Setup
-----------------------------

Have a :ref:`e3sm_diags development environment installed <dev-env>`
with the latest code from the master branch.


Running from the Cori Quick Start Guide
---------------------------------------

Below is the code to run :ref:`the parameter file linked here <cori-params>`,
call it ``lat_lon_demo.py``.

Remember to change ``prefix`` accordingly to a directory you can access.
Since we're running the software via Python, we can do Pythonic things, like using ``os.path.join()``.

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]

        prefix = '/global/project/projectdirs/acme/www/shaheen2/runs_with_api'
        param.results_dir = os.path.join(prefix, 'lat_lon_demo')

        runner.sets_to_run = ['lat_lon']
        runner.run_diags([param])

Run it like so with multiprocessing:

    ::

        python lat_lon_demo.py --multiprocessing --num_workers=32

We could have also set these multiprocessing parameters in the ``lat_lon_demo.py`` as well.
But we're showing that you can still submit parameters via the command line.

**Again, this new way of running is basically a replacement for the**
``params.py`` **in** ``e3sm_diags -p params.py -d diags.cfg`` **.**

Viewing the Results
^^^^^^^^^^^^^^^^^^^

`The results are here. <https://portal.nersc.gov/project/acme/shaheen2/dont_delete/lat_lon_demo/viewer/>`_
When you click on the Provenance folder and click on ``lat_lon_demo.py``, you can see a copy of the script used.

Again, in that folder, you can view the command used to run the software in ``cmd_used.txt``.
In this file, you can see that we ran with multiprocessing.

Also, look in the provenance under "Show Output Metadata".
Try modifying these values to change things.
You'll notice that to run the provenance, we have the set name (``lat_lon``) after the ``e3sm_diags`` argument.
This is a small change that a user shouldn't notice. Ex:

    ::

        e3sm_diags lat_lon --no_viewer --case_id 'CRU_IPCC' ...


Running with Custom Parameters for ``zonal_mean_2d``
----------------------------------------------------

In the file
`e3sm_diags/acme_diags/parameter/zonal_mean_2d_parameter.py <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/parameter/zonal_mean_2d_parameter.py>`_,
you can see the default value for ``plevs`` in the ``zonal_mean_2d`` set.
Say we want to change these to ``[10.0, 20.0, 30.0]``.

Call the file below ``zonal_mean_2d_plevs.py``. Notice that we are:

* Adding ``'zonal_mean_2d'`` to ``runner.sets_to_run``.
* Adding the new Parameter object, ``zonal_mean_2d_param``, to ``runner.run_diags()`` in addition to the original core parameter.
  * This is how the new ``plevs`` are propagated to the diags.

Also note that we can do more Pythonic things.
If we want to select only ``T`` and ``PRECT`` to run the diags on, you do what's commented out.
It's much easier than
:ref:`the other way, under "Using the selectors parameter" <selector-ex>`.

    .. code:: python

        import os
        from acme_diags.parameter.core_parameter import CoreParameter
        from acme_diags.parameter.zonal_mean_2d_parameter import ZonalMean2dParameter
        from acme_diags.run import runner

        param = CoreParameter()

        param.reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_e3sm_diags/climatology/'
        param.test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/climatology/'
        param.test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
        param.seasons = ["ANN"]

        prefix = '/global/project/projectdirs/acme/www/shaheen2/runs_with_api'
        param.results_dir = os.path.join(prefix, 'zonal_mean_2d_and_lat_lon_demo')

        # Uncomment the two lines below to just
        # run the diags with T and PRECT.
        # param.selectors += ['variables']
        # param.variables = ['T', 'PRECT']

        # The new changes are below.
        zonal_mean_2d_param = ZonalMean2dParameter()
        zonal_mean_2d_param.plevs = [10.0, 20.0, 30.0]

        runner.sets_to_run = ['zonal_mean_2d', 'lat_lon']
        runner.run_diags([param, zonal_mean_2d_param])


Run the diags:

    ::

        python zonal_mean_2d_plevs.py --multiprocessing --num_workers=32


Viewing the Results
^^^^^^^^^^^^^^^^^^^

`The results are located here. <https://portal.nersc.gov/project/acme/shaheen2/dont_delete/zonal_mean_2d_and_lat_lon_demo/viewer/>`_

Notice that though the ``plevs`` parameter was changed in the ``zonal_mean_2d`` plots,
it's unchanged for ``lat_lon``. For ``lat_lon``, you can see that we still have 200 and
850 as the pressure levels for certain diagnostics.
