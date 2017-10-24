
Installation, Configuration, and Running
========================================

Installation
------------

The installation procedure depends on what version you'd like to install.

Latest stable release
^^^^^^^^^^^^^^^^^^^^^

1. Make sure that you have Anaconda installed and upgraded to the latest version. If you don't have Anaconda
installed, look `here <https://conda.io/docs/user-guide/install/index.html#regular-installation>`_. 
Once installed, upgrade Anaconda like so:

   ::

       conda update conda

2. Get the yml file to create an environment. Use ``curl`` if on macOS.

   ::

       wget https://raw.githubusercontent.com/ACME-Climate/acme_diags/master/conda/acme_diags_env.yml

3. Remove any cached Anaconda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use Anaconda to create a new environment with ``acme_diags`` installed.  

   ::

       conda env create -f acme_diags_env.yml
       source activate acme_diags_env

   Tip: You can change the name of the environment using ``-n new_env_name``.


Environment for development
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Make sure Anaconda is installed and upgrade it to the latest version like so:

   ::

       conda update conda


2. Get the developmental yml file to create an environment.

   ::

       wget https://raw.githubusercontent.com/ACME-Climate/acme_diags/master/conda/acme_diags_env_dev.yml

3. Remove any cached Anaconda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use Anaconda to create a new environment. ``acme_diags`` **is not included in this environment.**

   ::

       conda env create -f acme_diags_env_dev.yml
       source activate acme_diags_env_dev

5. Get the latest code from master

   ::

       git clone https://github.com/ACME-Climate/acme_diags.git


   or if you already have a clone of the repo, pull the latest code from master.

   ::

       git pull origin master

5. Make and changes you want, then install.

   ::

       cd acme_diags
       python setup.py install

6. Run a quick test which generates one of each plot type. 
Remember to view the generated html which is here: ``all_sets/viewer/index.html``

   ::

       cd tests
       acme_diags_driver.py -d all_sets.cfg

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

