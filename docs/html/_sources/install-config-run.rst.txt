
Installation, Configuration, and Running
----------------------------------------

Installation
~~~~~~~~~~~~

Determining What Version to Get
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

-  **Latest stable release:**
-  If you have Anaconda and want to **create a new environment**:

   ::

       conda create -n acme_diags_env -c acme -c conda-forge -c uvcdat acme_diags

-  If you already have an environment:

   ::

       conda install -c acme -c conda-forge -c uvcdat acme_diags

-  **Nightlies: the latest code from master branch:** Go to the Anaconda
   `page for
   acme\_diags <https://anaconda.org/acme/acme_diags/files?channel=nightly>`__
   and choose a date, which is the version. For example, we'll choose
   June 23, 2017.
-  If you have Anaconda and want to **create a new environment**:

   ::

       conda create -n acme_diags_env -c acme/label/nightly -c conda-forge -c uvcdat acme_diags=2017.06.23

-  If you already have an environment:

   ::

       conda install -c acme/label/nightly -c conda-forge -c uvcdat acme_diags=2017.06.23

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

   or if you already have a clone of the repo

   ::

       git pull origin master

2. Make an environment using the Anaconda env file

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

Configuration
~~~~~~~~~~~~~

You must first do some configuration before you run the diagnostics.

1. Create a Python script, ex: ``myparams.py``. These scripts are simply
   just keys and values.
2. **At minimum, you must define values for the following:**

-  **``reference_data_path``: path to the reference (observational)
   data**
-  **``test_data_path``: path to the test (model output) data**
-  **``test_name``: name of the test (model output) file**

3. There are many other parameters that allow the user to customize
   regridding method, plotting backend, and much more. **A full list of
   parameters can be found `here <available-parameters.ipynb>`__.**

An example ``myparams.py`` script is shown below.

.. code:: ipython2

    reference_data_path = '/Users/zhang40/Documents/AIMS/amwg/amwg20140804/obs_data_20140804/'
    test_data_path = '/Users/zhang40/Documents/ACME_simulations/'
    
    test_name = '20160520.A_WCYCL1850.ne30'
    
    # a few optional parameters
    backend = 'vcs'  # use 'mpl' for matplotlib
    
    results_dir = 'myresults'  # name of folder where all results will be stored

To add your own diagnostics, create a cfg file like the one below. **All
of the keys in the cfg file are possible parameters as well. A full list
is `here <available-parameters.ipynb>`__**. We'll call this
``mydiags.cfg``.

.. code:: ipython2

    ```
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
    ```

If you have multiple diagnostics you want to run, create a cfg file like
the one below.

::

    [Diags]
    # put all of the parameters for a diags run here

    [Diags 2]
    # another diags run
    # make sure that the title ("Diags 2") is unique.

Running
~~~~~~~

If you **don't** have your own diagnostics, simply just run:

``acme_diags_driver.py -p myparams.py``

If you do have your own own diagnostics, run:

``acme_diags_driver.py -p myparams.py -d mydiags.cfg``

View the results by opening ``index.html`` in the location specified.
