Defining Parameters
===================

Ways to define parameters
-------------------------

There are three ways to input parameters to the diagnostics: 

1. **Command line**: For example: ``e3sm_diags -p myparam.py --variables T PRECT`` 
   will set the variables to ``['T', 'PRECT']``. 
2. **Parameters file**: In the command ``e3sm_diags -p myparam.py``, 
   the parameters file is ``myparam.py``. 
3. **Diagnostics file**: In the command ``e3sm_diags -d mydiags.cfg``, 
   the diagnostics file is ``mydiags.cfg``.

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

    [#]
    variables = ["PRECT"]
    regions = ["global"]
    seasons = ["ANN"]

    [#]
    variables = ["SST"]
    regions = ["ocean"]
    seasons = ["JJA"]

Running with just the parameters file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Running ``e3sm_diags -p myparams.py`` will just run the
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

Running ``e3sm_diags -d mydiags.cfg`` will have two runs with
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

Running ``e3sm_diags -p myparams.py -d mydiags.cfg`` will also
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
``e3sm_diags -p myparams.py -d mydiags.cfg --variables PREH2O``
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

Selecting certain parameters
----------------------------

When you run ``e3sm_diags`` with a file passed in via ``-p``,
the parameters in that file are inserted into each diagnostics from files like
`this <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/default_diags/lat_lon_model_vs_obs.cfg/>`_,
overwriting any duplicates in the process.
A single diagnostics starts with ``[#]``.
If you provide your own cfg file with ``-d``, the same happens.

For example, say we have the following parameters in a Python file:

.. code:: python

    reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
    test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    variables = ['PRECT']
    sets = ['lat_lon']

Since we're running the ``lat_lon`` plotset, and since it defaults to ``model_vs_obs``,
it will open 
`this file <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/default_diags/lat_lon_model_vs_obs.cfg/>`_.
Each of the parameters in the Python file will be inserted into each of the diagnostics runs.


So each of the 100+ ``lat_lon`` diagnostics will be done with ``variables = ['PRECT']``.
However, this is nonsensical.
**What we want to do is to "select" the diagnostics** `from here <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/driver/default_diags/lat_lon_model_vs_obs.cfg/>`_ **that use PRECT.**

Using the selectors parameter
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the above Python file, we can designate the ``variables`` parameter to be a "selector".

First, find the default ``selectors`` used
`here <https://github.com/E3SM-Project/e3sm_diags/blob/master/acme_diags/acme_parameter.py>`_
and copy what current parameters are used as selectors. It's what's defined by ``self.selectors``.

In **your Python file**, paste these along with any parameters you want as selectors.
It should look something like this:

.. code:: python

    reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
    test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    variables = ['PRECT']
    sets = ['lat_lon']
    # 'sets' and 'seasons' were our default values.
    selectors = ['sets', 'seasons', 'variables']


If you run ``e3sm_diags`` now like this, you'll only run the diagnostics
that had the variables originally as ``'PRECT'``.

**Remember that you can use any of the parameters defined**
`here <available-parameters.html>`_
**as selectors.**

Say we only wanted to select the diagnostics using ``'PRECT'`` and using specific observational data. 
We can do the following:

.. code:: python

    reference_data_path = '/global/project/projectdirs/acme/acme_diags/obs_for_acme_diags/'
    test_data_path = '/global/project/projectdirs/acme/acme_diags/test_model_data_for_acme_diags/'

    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'

    variables = ['PRECT']
    sets = ['lat_lon']
    ref_name = ['GPCP_v2.2', 'ERA-Interim']
    # 'sets' and 'seasons' were our default values.
    selectors = ['sets', 'seasons', 'variables', 'ref_name']
