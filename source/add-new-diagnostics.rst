
How to Add New Diagnostics
==========================

The following guide will explain how to add new diagnostics. There is
also an example, detailing a sample process from start to finish.

Adding New Diagnostics
----------------------

**Before you do this, please *actually* read `Defining Parameters and
All Available Parameters <available-parameters.ipynb>`__. We're using
json files to add multiple diagnostics.**

Examples of parameters
~~~~~~~~~~~~~~~~~~~~~~

The following examples are of common parameters used in diagnostics.
Again, for more parameters and detail, please refer to the Defining
Parameters and All Available Parameters linked previously.

When creating json files, it's recommended that you check your syntax
using a linter (like `jsonlint <http://jsonlint.com/>`__) to make sure
that it's valid. Or better yet, just use .cfg files.

-  **The name of the folder where all of the output is created.**

   ::

       results_dir= "myresults"

-  **Each run must have a ``case_id``, which is the name of the folder
   where the outputs will be stored relative to the path of
   ``results_dir``. So for this run, the directory structure is:
   ``myresults/set5_MERRA``.**

   ::

       case_id = "set5_MERRA"

-  **The ``sets`` is a list of sets ran. The actual numbers can be
   integers or strings. You can choose one or more. We currently support
   AMWG sets 3, 4, 5, 7, 13. Equivalently values are 'zonal\_mean\_xy',
   'zonal\_mean\_2d', 'lat\_lon', 'polar', and 'cosp\_histogram'.**

   ::

       sets = ['lat_lon']

-  **Add in the name of the observations you're going to use.**

   ::

       ref_name = "MERRA"

-  **A list of all variables to use.**

   ::

       variables = ["T"]

-  **Regions are added like so. If no value for ``regions`` is given,
   it's ``global`` by default. A full list of all regions is
   `here <https://github.com/ACME-Climate/acme_diags/blob/master/acme_diags/derivations/default_regions.py>`__**

   ::

       regions = ["land", "ocean_TROPICS"]

-  **The following seasons are supported.**

   ::

       seasons = ["ANN", "DJF", "MAM", "JJA", "SON"]

-  **For 3-d variables, pressure levels can be added.**

   ::

       plevs = [200.0, 850.0]

Adding derived variables
~~~~~~~~~~~~~~~~~~~~~~~~

We have a set of built-in derived variables for the ACME model
diagnostics
`here <https://github.com/ACME-Climate/acme_diags/blob/master/acme_diags/derivations/acme.py>`__
(search for ``derived_variables``). The diagnostics software looks into
the ``derived_variables`` dictionary for variable keys and operations
needed for deriving new variables (renaming, unit conversions,
calculations, etc).

If users want to, they can add their own derived variables, which is
added to the default list during runtime and overwrite any default
values if there's a collision. Since derived variables require code,
such functionality cannot be added to json/cfg files. We do the
following in the parameters script, which is a Python script (ex: the
Python script is ``myparams.py`` in the command
``acme_diags_driver.py -p myparams.py -d mydiags.cfg``).

Format of the ``derived_variables`` dictionary
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code:: python

    derived_variables = {
        'user_inputted_var': {
            ('output_var1', 'output_var2'): function_to_call
        }
    }

Above is how a ``derived_variables`` dictionary is formatted.
``'user_inputted_var'`` is the variable that is defined in the
``variables`` part of the json file. ``'output_var1'`` and
``'output_var2'`` are the variables inside the ``test_name`` (model)
file.

**Example of adding derived variables to a parameters script**

.. code:: python

    # myparams.py

    def albedo_obs(rsdt, rsut):
        """ TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension """
        var = rsut / rsdt
        var.units = "dimensionless"
        var.long_name = "TOA albedo"
        return var

    derived_variables = {
        'ALBEDO': {
            ('rsdt', 'rsut'): albedo_obs
        }
    }

The above code will allow it so that if ``variables = "ALBEDO"`` in the
diagnostics file (the json or cfg file) and ``rsdt`` and ``rsut`` are
variables in the test (model) file, the ``albedo_obs()`` function is ran
on the ``rsdt`` and ``rsut`` variables from the test (model) file.

Example
-------

The example below will do one diagnostics run globally with the
``ALBEDO`` variable, annually. Below is the json file, call it
``mydiags.cfg``.

::

    [Diags]
    sets = ['lat_lon']
    case_id = "lat_lon_CERES"
    variables = ["ALBEDO"]
    ref_name = "edition_4_ceres_ebaf_toa"
    reference_name = "edition_4_ceres_ebaf_toa"
    seasons = ["ANN"]
    regions = ["global"]
    contour_levels = [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95]
    diff_levels = [-0.25, -0.2, -0.15, -0.1, -0.07, -0.05, -0.03, 0.0, 0.03, 0.05, 0.07, 0.1, 0.15, 0.2, 0.25]

And below is the parameters file, call it ``myparams.py``. This is
tested to run on aims4. To run on another machine, please edit the
``reference_data_path``, ``test_data_path``, and ``test_name``
accordingly.

.. code:: python

    reference_data_path = '/space1/test_data/CERES-EBAF/'
    test_data_path = '/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

    test_name = '20160520.A_WCYCL1850.ne30'

    backend = 'vcs'
    diff_title = 'Test - Reference'
    results_dir = 'myresults'

    def albedo_obs(rsdt, rsut):
        """TOA (top-of-atmosphere) albedo, (solin - fsntoa) / solin, unit is nondimension"""
        var = rsut / rsdt
        var.units = "dimensionless"
        var.long_name = "TOA albedo"
        return var

    derived_variables = {
        'ALBEDO': {
            ('rsdt', 'rsut'): albedo_obs
        }
    }

Run the command like so:
``acme_diags_driver.py -p myparams.py -d mydiags.cfg``
