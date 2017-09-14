
Installation, Configuration, and Running
========================================

Installation
------------

The installation procedure depends on what version you'd like to install.

Latest stable release
^^^^^^^^^^^^^^^^^^^^^

-  If you have Anaconda and want to **create a new environment**:

   ::

       conda create -n acme_diags_env -c acme -c conda-forge -c uvcdat acme_diags

-  If you want to install in an existing environment:

   ::

       conda install -c acme -c conda-forge -c uvcdat acme_diags


Nightlies: the latest code from master branch
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Go to the Anaconda 
`page for acme\_diags <https://anaconda.org/acme/acme_diags/files?channel=nightly>`__
and choose a date, which is the version. For example, we'll choose
Sept 14, 2017.

-  If you have Anaconda and want to **create a new environment**:

   ::

       conda create -n acme_diags_env -c acme/label/nightly -c conda-forge -c uvcdat acme_diags=2017.09.14

-  If you want to install in an existing environment:

   ::

       conda install -c acme/label/nightly -c conda-forge -c uvcdat acme_diags=2017.09.14

-  **Optional for vcs:** If you plan on using ``vcs`` and you don't want
   to use the X windowing system, run the following command:

   ::

       conda install mesalib -c conda-forge -c uvcdat


Environment for development
^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you want to develop/modify the code for ``acme_diags``, follow these
steps:

1. Get the latest code of acme\_diags

   ::

       git clone https://github.com/ACME-Climate/acme_diags

   or if you already have a clone of the repo, pull the latest version:

   ::

       git pull origin master

2. Create an environment using the Anaconda env file

   ::

       conda env create -f acme_diags/conda/acme_diags_env.yml

   Tip: if you want to change the name of the env, just append the
   following: ``-n new_name``

3. Activate the environment using whatever name you used

   ::

       source activate acme_diags_env

4. Proceed to make changes to your code, then go to the location of
   ``setup.py`` and do the following

   ::

       python setup.py install

   Note that we you'll need to repeat the installation step every time you make
   changes to the source code and want them to take effect.

Configuration
-------------

You must first do some configuration before you run the diagnostics.

1. Create a Python script, ex: ``myparams.py``. This script contains simply
   pairs of keys and values. **At minimum, you must define values for the following:**

-  **``reference_data_path``: path to the reference (observational)
   data**
-  **``test_data_path``: path to the test (model output) data**
-  **``test_name``: name of the test (model output) file**

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

    [Diags]
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

    [Diags]
    # put all of the parameters for a diags run here

    [Diags 2]
    # another diags run
    # make sure that the title ("Diags 2") is unique.


Running
-------

If you **don't** have your own diagnostic file (e.g. ``mydiags.cfg``), simply run: ::

  acme_diags_driver.py -p myparams.py

to generate the standard set of ACME diagnostics figures.
If you do have your own own diagnostic file, specify it on the command line: ::

  acme_diags_driver.py -p myparams.py -d mydiags.cfg

View the results by opening ``index.html`` in the location specified.

