Configuration and Running v1
============================
This guide is made for ``e3sm_diags`` v1. We are preparing for a new document for v2. Although this guide still works, for the time being if you are using ``e3sm_diags`` v1, we recommend you to use the quick guides to start a run.

Configuration
-------------

You must first do some configuration before you run the diagnostics.

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

    reference_data_path = '/p/cscratch/acme/data/obs_for_acme_diags/'
    test_data_path = '/p/cscratch/acme/data/test_model_data_for_acme_diags/'
    
    test_name = '20161118.beta0.FC5COSP.ne30_ne30.edison'
    
    # a few optional parameters
    backend = 'vcs'  # use 'mpl' for matplotlib
    
    results_dir = 'myresults'  # name of folder where all results will be stored

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
