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
