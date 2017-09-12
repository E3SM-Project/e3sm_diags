
Defining Parameters and All Available Parameters
================================================

Ways to define parameters
-------------------------

There are three ways to input parameters to the diagnostics: 1. Command
line: For example:
``acme_diags_driver.py -p myparam.py --variables T PRECT`` will set the
variables to ``['T', 'PRECT']``. 2. Parameters file: In the command
``acme_diags_driver.py -p myparam.py``, the parameters file is
``myparam.py``. 3. Diagnostics file: In the command
``acme_diags_driver.py -d mydiags.cfg``, the diagnostics file is
``mydiags.cfg``.

**Each of these ways have a level of priority, with the command line
input having the highest priority and the diagnostics file having the
lowest priority.**

Examples
~~~~~~~~

Say that we have the following files:

``myparams.py:``

.. code:: python

    sets = ['lat_lon']
    variables = ['T']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

``mydiags.cfg:``

::

    [Diags]
    variables = ["PRECT"]
    regions = ["global"]
    seasons = ["ANN"]

    [Diags 2]
    variables = ["SST"]
    regions = ["ocean"]
    seasons = ["JJA"]

Running with just the parameters file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running ``acme_diags_driver.py -p myparams.py`` will just run the
lat-lon contour diagnostics once with the parameters being:

.. code:: python

    variables = ['T']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

Running with the diagnostics file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The ``cfg`` files support all of the non-programatic parameters covered
below in the Available Parameters section. ``cfg`` files are also the
way of defining multiple diagnostics runs, in ``mydiags.cfg`` defined
above, we have two runs.

Running ``acme_diags_driver.py -d mydiags.cfg`` will have two runs with
the following parameters

Run 1:

.. code:: python

    variables = ['PRECT']
    regions = ['global']
    seasons = ['ANN']

Run 2:

.. code:: python

    variables = ['SST']
    regions = ['ocean']
    seasons = ['JJA']

Running ``acme_diags_driver.py -p myparams.py -d mydiags.cfg`` will also
have two runs, but the parameters in ``myparams.py`` will take priority
over the ones described in ``mydiags.cfg``. So the runs will be:

Run 1:

.. code:: python

    variables = ['T']
    regions = ['global']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

Run 2:

.. code:: python

    variables = ['T']
    regions = ['ocean']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

Running with command line arguments
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

With all of the three aforementioned ways of setting the parameters
(just ``myparams.py``, just ``mydiags.cfg``, or both ``myparams.py`` and
``mydiags.cfg``), command line arguments can be added to all.

So running
``acme_diags_driver.py -p myparams.py -d mydiags.cfg --variables PREH2O``
will have the variables in both runs be ``PREH2O``:

Run 1:

.. code:: python

    variables = ['PREH2O']
    regions = ['global']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

Run 2:

.. code:: python

    variables = ['PREH2O']
    regions = ['ocean']
    seasons = ['DJF', 'MAM', 'JJA', 'SON']

--------------

Available Parameters
--------------------

The driver needs a parameters file to run. In these files, there is
support for many features related to diagnostics.

Given a command like ``acme_diags_driver.py -p params.py``, the
parameters in ``params.py`` will overwrite any predefined values for all
of the runs.

Parameters for diagnostics
~~~~~~~~~~~~~~~~~~~~~~~~~~

Any parameter that is red is one that is required for a diagnostics run.
The rest are optional and often have default values.

The parameters in the section are ones related to test/reference
specifications. Below are the parameters related to file I/O.

-  **``results_dir``**: the name of the folder where all runs will be
   stored. If not defined, the folder where all of the results are
   created in is named ``acme_diags_results-TIMESTAMP``.
-   **``case_id``**: the name of the folder where the results (nc files
   and plots) will be stored for a single run. ex:
   ``results_dir/case_id``\ 
-   **``reference_data_path``**: path to the reference (obs) data.
-   **``test_data_path``**: path to the test (model) data.
-  **``reference_name``**: the name of the reference (obs) file. This
   doesn't need to be defined if your running the default AMWG sets. In
   the built-in parameters files for these, the ``reference_name`` is
   already defined.
-  \ **``test_name``**: the name of the test (model output) file.
-  **``save_netcdf``**: set to ``True`` if you want the reference, test,
   and difference data saved. It's ``False`` by default.

The parameters below are for running the diagnostics in parallel using
multiprocessing or distributedly.

-  **``num_workers``**: number of processes to use. If not defined, it
   will used the number of logical (not physical) processes on the
   machine. ex: ``num_workers = 8``
-  **``multiprocessing``**: set to ``True`` to use multiprocessing. It's
   ``False`` by default. ``multiprocessing`` and ``distributed`` cannot
   both be ``True.``
-  **``distributed``**: set to ``True`` to run the diagnostics
   distributedly. It's ``False`` by default. ``multiprocessing`` and
   ``distributed`` cannot both be ``True.``

The parameters below are related to the actual climate-related
functionality of the diagnostics.

-  \ **``sets``**: A list of the sets to be ran. Ex.
   ``sets=['zonal_mean_xy', 'zonal_mean_2d', 'lat_lon, 'polar', 'cosp_histogram']``
   or ``sets=['3', '4', '5, '7', '13']``\ 
-  **``datasets``**: A list of what datasets to use. Based on this and
   ``sets``, it will load the corresponding ``*ACME.json`` files located
   `here <https://github.com/ACME-Climate/acme_diags/tree/master/acme_diags/driver>`__.
   Possible values are: "ACME" or "AMWG". It's ``datasets=['ACME']`` by
   default when nothing is defined.
-  \ **``variables``**: What variables to use for this run. Ex:
   ``variables=["T", "PRECT"]``.
-  \ **``seasons``**: A list of season to use. Possible values are:
   "ANN", "DJF", "MAM", "JJA", "SON". Ex:
   ``seasons=["ANN", "DJF", "MAM", "JJA", "SON"]``.
-  **``regions``**: A list of regions. If none, it's ``global`` by
   default. See `this
   file <https://github.com/zshaheen/acme_diags/blob/master/acme_diags/derivations/default_regions.py>`__
   for a list of possible regions. Ex: ``regions=["global","TROPICS"]``.
-  **``plevs``**: A list of pressure levels to use. Ex:
   ``plevs=[850.0, 200.0]``.
-  **``regrid_tool``**: The regrid tool to use. **Possible values are:
   \_\_\_\_.** Is ``'esmf'`` by default when no value is given.
-  **``regrid_method``**: What regird method of the regrid tool to use.
   **Possible values are ``'linear'``, ``'conservative'``.** Is
   ``'linear'`` by default when no value is given.

Parameters for plotting
~~~~~~~~~~~~~~~~~~~~~~~

The figure below is an sample output. We use this to described what each
plotting parameter does.

.. figure:: parameter_example.png
   :alt: Example

   Figure1

Below are general plotting-related parameters.

-  **``main_title``**: Main title of the image. It's "PRECT ANN global"
   in the example and is blank by default.
-  **``backend``**: Can either be ``vcs`` or
   ``cartopy``/``mpl``/``matplotlib``. Is ``vcs`` by default when no
   value is given.
-  **``output_file``**: Name of the output file. Is ``output`` by
   default. Remember **not** to include an extension.
-  **``output_format``**: A list of formats that yout want the plot to
   be output to. Can be something like ``['png', 'pdf', 'svg'].``\ Is
   ``['png']`` when nothing is present.
-  **``canvas_size_w``**: width of the image in pixels and only used by
   vcs. Is 1212 by default.
-  **``canvas_size_h``**: height of the image in pixels and only used by
   vcs. Is 1628 by default.
-  **``arrows``**: Is either ``True`` (default value) or ``False`` and
   will accordingly show or hide the arrows on the legend for all of the
   graphs.
-  **``logo``**: ``True`` (default value) to show the UV-CDAT logo on
   the vcs backend, ``False`` to not. Just keep it on please.

The parameters below are for each of the three plots (``test``,
``reference``, and ``diff``) in the image.

-  **``test_title``**: the title for the test plot. It's "Test Title" in
   the image and is blank by default.
-  **``test_colormap``**: If not defined in the parameters, the default
   value is ``'cet_rainbow.rgb'``. Matplotlib colormaps are supported.
   Users can even use colormaps located in
   ``acme_diags/plot/colormaps/``, by referencing them by the filename
   (ex: ``'cet_rainbow.rgb'``). Also, paths to a custom ``.rgb`` file is
   supported.
-  **``contour_levels``**: the levels on the legend of the test and
   reference plot. It's [0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12,
   14, 17] in the image but automatically gets the range by default.
-  **``test_units``**: the units that are on the top-right of the test
   plot. It's "mm/day" in the image. If blank, it automatically gets the
   units from the test data.

The ``reference`` and ``diff`` plots also have the same keywords which
are semantically the same for their respective plots. Below are the
values they hold for the image above.

-  **``reference_title``**: "Reference title" in the image and is blank
   by default.
-  **``reference_colormap``**: This wasn't define in the parameter, so
   it is the default value, ``'cet_rainbow.rgb'``. Matplotlib colormaps
   are supported. Users can even use colormaps located in
   ``acme_diags/plot/colormaps/``, by referencing them by the filename
   (ex: ``'cet_rainbow.rgb'``). Also, paths to a custom ``.rgb`` file is
   supported.
-  **``contour_levels``**: You only need one ``contour_levels`` in you
   script. It's used in the reference plot. [0, 0.2, 0.5, 1, 2, 3, 4, 5,
   6, 7, 8, 9, 10, 12, 14, 17] in the image.
-  **``reference_units``**: "mm/day" in the image. If blank, it
   automatically gets the units from the reference data.

-  **``diff_title``**: "Test - Reference" in the image. If blank, the
   default is "Model - Observation".
-  **``diff_colormap``**: is ``'RdBu_r'`` in the image above and
   ``'bwr'`` by default. Matplotlib colormaps are supported. Users can
   even use colormaps located in ``acme_diags/plot/colormaps/``, by
   referencing them by the filename (ex: ``'cet_rainbow.rgb'``). Also,
   paths to a custom ``.rgb`` file is supported.
-  **``diff_levels``**: [-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3,
   4, 5, 6] in the image but automatically gets the range by default.
-  **``diff_units``**: "mm/day" in the image. If blank, it automatically
   gets the units from the test - reference data.
