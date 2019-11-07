
Installation
============

The installation procedure depends on what version you'd like to install.

Activate **e3sm_unified** environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have an account on one of the E3SM supported machines (Cori, Compy, Acme1, Anvil, Cooley, Rhea), you
can access ``e3sm_diags`` by activating ``e3sm_unified``, which is a conda environment pulls together Python
and other E3SM analysis tools such as ``e3sm_diags``, ``mpas-analysis``, ``NCO``, ``cdat`` and ``processflow``.

The paths to ``e3sm_unified`` activation scripts are machine depended:

**Compy**
    ::

     source /compyfs/software/e3sm-unified/load_latest_e3sm_unified.sh


**Cori**
    ::

     source /global/project/projectdirs/acme/software/anaconda_envs/load_latest_e3sm_unified.sh
    

**Anvil/blues**
    ::

     source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified.sh


**Cooley**
    ::

     source /lus/theta-fs0/projects/ccsm/acme/tools/e3sm-unified/load_latest_e3sm_unified.sh


**acme1**
    ::

     source /usr/local/e3sm_unified/envs/load_latest_e3sm_unified.sh


**Rhea**
    ::

     source /ccs/proj/cli900/sw/rhea/e3sm-unified/load_latest_e3sm_unified.sh
 

Change ``.sh`` to ``.csh`` for csh shells.
Note that ``e3sm_unified``'s development cycle is not in phase with ``e3sm_diags``, therefore the version of ``e3sm_diags`` included may not be the latest. To install latest stable releases, refer to following:

.. _install_latest:

Latest stable release
^^^^^^^^^^^^^^^^^^^^^

1. Make sure that you have Anaconda installed and upgraded to the latest version. If you don't have Anaconda
installed, look `here <https://docs.conda.io/projects/conda/en/latest/user-guide/install/>`_. 
Once installed, upgrade Anaconda like so:

   ::

       conda config --set ssl_verify false
       binstar config --set ssl_verify False
       conda update conda

2. Get the yml file to create an environment. Use ``curl`` if on macOS.

   ::

       wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env.yml

3. Remove any cached Anaconda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use Anaconda to create a new environment with ``e3sm_diags`` installed.  

   Tip: You can change the name of the environment by adding ``-n new_env_name`` to the end of `conda env create ...`.

   ::

       conda env create -f e3sm_diags_env.yml
       source activate e3sm_diags_env


.. _dev-env:

Environment for development
^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Make sure Anaconda is installed and upgrade it to the latest version like so:

   ::

       conda config --set ssl_verify false
       binstar config --set ssl_verify False
       conda update conda


2. Get the developmental yml file to create an environment.

   ::

       wget https://raw.githubusercontent.com/E3SM-Project/e3sm_diags/master/conda/e3sm_diags_env_dev.yml

3. Remove any cached Anaconda packages. This will ensure that you always get the latest packages.

   ::

       conda clean --all

4. Use Anaconda to create a new environment. ``e3sm_diags`` **is not included in this environment.**

   ::

       conda env create -f e3sm_diags_env_dev.yml
       source activate e3sm_diags_env_dev

5. Get the latest code from master

   ::

       git clone https://github.com/E3SM-Project/e3sm_diags.git


   or if you already have a clone of the repo, pull the latest code from master.

   ::

       git pull origin master

5. Make and changes you want, then install.

   ::

       pip install .

6. Run a quick test which generates one of each plot type. 
Remember to view the generated html located here: ``all_sets/viewer/index.html``

   ::

       cd tests/system
       python all_sets.py -d all_sets.cfg
