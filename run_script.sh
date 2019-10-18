# module load python/2.7-anaconda-4.4
# source activate e3sm_diags_env_dev
pip install /global/homes/f/forsyth/e3sm_diags/

rm -r /global/project/projectdirs/acme/www/forsyth/enso_diags_1_8
python run_enso_diags.py -d enso.cfg
chmod -R 755 /global/project/projectdirs/acme/www/forsyth/enso_diags_1_8
echo 'Visit https://portal.nersc.gov/project/acme/forsyth/enso_diags_1_8'
