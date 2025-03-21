
General quick guide for running e3sm_diags v3
=========================================================================


1. Installation
-----------------------------------------------------------

We will use the e3sm_unifed environment to install.
For the latest stable release or if you don't have access to e3sm analysis machines,
please instead refer to :ref:`Latest stable release <install_latest>`.

Most of the E3SM analysis software is maintained with an Anaconda metapackage
(E3SM unified environment).
If you have an account on an E3SM supported machine
(**Andes, Anvil, Chrysalis, Compy, Perlmutter**),
then to get all of the tools in the metapackage in your path,
use the corresponding activation command below.
(Change ``.sh`` to ``.csh`` for csh shells.)

Below, we also provide the paths for observational data needed by ``e3sm_diags`` (<obs_path>),
and some sample model data for testing (<test_data_path>).
Both <obs_path> and <test_data_path> have two subdirectories:
``/climatology`` and ``/time-series`` for climatology and time-series data respectively.

Also listed below are paths where the HTML files (<html_path>) must be located to be displayed
at their corresponding web addresses (<web_address>).
Note that only some machines (**Acme1, Anvil, Chrysalis, Compy, Cori**) have a web server.


Andes
^^^^^
<activation_path>: ``source /ccs/proj/cli900/sw/rhea/e3sm-unified/load_latest_e3sm_unified_andes.sh``

<obs_path>:``/ccs/proj/cli115/e3sm_diags_data/obs_for_e3sm_diags/``

<test_data_path>:``/ccs/proj/cli115/e3sm_diags_data/test_model_data_for_acme_diags/``

Anvil
^^^^^
<activation_path>: ``source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_anvil.sh``

<obs_path>: ``/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/``

<test_data_path>: ``/lcrc/soft/climate/e3sm_diags_data/test_model_data_for_acme_diags/``

<html_path>: ``/lcrc/group/e3sm/public_html/diagnostic_output/<username>/``

<web_address>: ``https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/<username>``

Chrysalis
^^^^^^^^^
<activation_path>: ``source /lcrc/soft/climate/e3sm-unified/load_latest_e3sm_unified_chrysalis.sh``

<obs_path>: ``/lcrc/soft/climate/e3sm_diags_data/obs_for_e3sm_diags/``

<test_data_path>: ``/lcrc/soft/climate/e3sm_diags_data/test_model_data_for_acme_diags/``

<html_path>: ``/lcrc/group/e3sm/public_html/diagnostic_output/<username>/``

<web_address>: ``https://web.lcrc.anl.gov/public/e3sm/diagnostic_output/<username>``


Compy
^^^^^
<activation_path>: ``source /share/apps/E3SM/conda_envs/load_latest_e3sm_unified_compy.sh``

<obs_path>: ``/compyfs/e3sm_diags_data/obs_for_e3sm_diags/``

<test_data_path>: ``/compyfs/e3sm_diags_data/test_model_data_for_acme_diags/``

<html_path>: ``/compyfs/www/<username>``

<web_address>: ``https://compy-dtn.pnl.gov/<username>``



Perlmuter
^^^^^^^^^^^^
<activation_path>: ``source /global/cfs/cdirs/e3sm/software/anaconda_envs/load_latest_e3sm_unified_pm-cpu.sh``

<obs_path>: ``/global/cfs/cdirs/e3sm/e3sm_diags/obs_for_e3sm_diags/``

<test_data_path>: ``/global/cfs/cdirs/e3sm/e3sm_diags/test_model_data_for_acme_diags/``

<html_path>: ``/global/cfs/cdirs/e3sm/www/<username>``

<web_address>: ``http://portal.nersc.gov/cfs/e3sm/<username>``

2. Config and run
--------------------------------------------------------

See :doc:`Anvil <quick-guide-anvil>`, `Chrysalis <quick-guide-chrysalis>`, :doc:`Compy <quick-guide-compy>`,
or :doc:`Perlmutter <quick-guide-perlmutter>`. If you are using other machines, like Andes,
then follow one of these guides, substituting the corresponding paths from above.
Again, note that only some machines (**Anvil, Compy, Chrysalis and Perlmutter**) have a web server.

