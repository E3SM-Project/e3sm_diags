Installation
============

The installation procedure depends on what version you'd like to install.

Activate **e3sm_unified** environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
If you have an account on one of the E3SM supported machines (NERSC, Compy, LCRC, Cooley, Rhea), you
can access ``e3sm_diags`` by activating ``e3sm_unified``, which is a conda environment that pulls together Python
and other E3SM analysis tools such as ``e3sm_diags``, ``mpas-analysis``, ``NCO``, and ``processflow``.

The paths to ``e3sm_unified`` activation scripts are machine dependent:

**Compy**
    ::

     source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh


**NERSC Perlmutter**
    ::

     source /global/common/software/e3sm/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh

**LCRC**
    ::

     source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh


**Andes**
    ::

     source /ccs/proj/cli115/software/e3sm-unified/load_latest_e3sm_unified_andes.sh

Change ``.sh`` to ``.csh`` for csh shells.
Note that ``e3sm_unified``'s development cycle is not in phase with ``e3sm_diags``,
therefore the version of ``e3sm_diags`` included may not be the latest.
To install latest stable releases, refer to following:

.. _conda_environment:

Installation in a Conda Environment
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If the E3SM Unified environment doesn't serve your needs, you can alternatively
install the latest version in your own custom conda environment.

First, activate conda or install it if not available. Examples shown below, details vary on the machine.

Compy
~~~~~
    ::

     module load anaconda3/2019.03
     source /share/apps/anaconda3/2019.03/etc/profile.d/conda.sh


NERSC
~~~~~
    ::

     module load python/3.7-anaconda-2019.10
     source /global/common/cori_cle7/software/python/3.7-anaconda-2019.10/etc/profile.d/conda.sh

.. _conda_environment_others:

Others/Local
~~~~~~~~~~~~

In most cases, a user-installed conda is preferred, especially if you need an up-to-date, isolated Conda environment without affecting system packages., follow these instructions for Unix-like platforms (macOS & Linux)

1. Download Miniforge

    ::

        wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"

2. Install Miniforge

    ::

        bash Miniforge3-$(uname)-$(uname -m).sh

    When you see: ::

        by running conda init? [yes|no]
        [no] >>> yes

    respond with ``yes`` so ``conda`` commands are available on
    initializing a new bash terminal.

3. If you are working on a machine/network that intercepts SSL communications, you will get
an SSL error unless you disable the SSL verification:

    ::

        conda config --set ssl_verify false
        binstar config --set ssl_verify False


4. Once conda are properly working, you can install the **(a) Latest Stable Release** or
create a **(b) Development Environment**.

.. _install_latest:

(a) Latest Stable Release
-------------------------

1. Follow :ref:`"Others/Local" <conda_environment_others>` section for installing conda.

Create a new conda environment with ``e3sm_diags`` installed and activate it:

    ::

        # Tip: Add the flag ``-n <name_of_env>`` to customize the name of the environment
        conda create -n e3sm_diags_env e3sm_diags
        conda activate e3sm_diags_env

.. _dev-env:

(b) Development Environment
---------------------------

1. Clone e3sm_diags:

    ::

    	git clone https://github.com/E3SM-Project/e3sm_diags.git

2. Enter the e3sm_diags directory.

    ::

        cd e3sm_diags

3. Use conda to create a new dev environment (``e3sm_diags`` **is not included in this environment**).

    - Tip: Add the flag ``-n <name_of_env>`` to customize the name of the environment

    ::

        conda env create -f conda-env/dev.yml
        conda activate e3sm_diags_env_dev

4. Unlike the latest stable release (i.e., the user environment), the development environment does not include E3SM Diags (``e3sm-diags``). Instead, the developer will ``make install`` (or ``python -m pip install .``) to build ``e3sm-diags`` with changes. Note: Before rebuilding, always run ``make clean`` to remove old build files. This ensures a correct and clean reinstallation. 


.. _dev-env-long:

(c) Development Environment (Longer version for contributing code)
------------------------------------------------------------------

.. note::
    The dev environment includes quality assurance (QA) tools such as code formatters,
    linters, and ``pre-commit``. **You will need to use the dev environment for all
    contributions** because these QA tools are enforced using ``pre-commit`` checks in
    the continuous integration/continuous deployment build.

1. Follow :ref:`"Others/Local" <conda_environment_others>` section for installing conda.

2. Create a new fork of e3sm_diags:

    ::

        # Go to https://github.com/E3SM-Project/e3sm_diags
        # Click "Fork" in the upper right hand corner. This will fork the main repo.
        # Click the green "Code" button
        # Choose the HTTPS or SSH option.
        # Click the clipboard icon to copy the path.
        # On your command line:
        git clone <path>
        git remote -v
        # You should see your fork listed as `origin`



    (Optional) add the main e3sm_diags repository as an upstream remote:
    ::

        git remote add upstream https://github.com/E3SM-Project/e3sm_diags.git
        # You're now ready to start working on your fork.
     

3. Remove any cached conda packages. This will ensure that you always get the latest packages.

    ::

        conda clean --all

4. Enter the fork directory.

    ::

        cd e3sm_diags

5. Use conda to create a new dev environment (``e3sm_diags`` **is not included in this environment**).

    - Tip: Add the flag ``-n <name_of_env>`` to customize the name of the environment

    ::

        conda env create -f conda-env/dev.yml
        conda activate e3sm_diags_env_dev

6. Install ``pre-commit``.

    ::

        pre-commit install

7. Make the desired changes to E3SM Diags, then rebuild and install with:

    ::

        make install # or python -m pip install .

8. Check that tests pass: ``./tests/test.sh``. This takes about 4 minutes.

9. Commit changes and make sure ``pre-commit`` checks pass
    ::

        git commit -m "..."

    .. figure:: pre-commit-passing.png
       :alt: pre-commit Output

       ``pre-commit`` Output
