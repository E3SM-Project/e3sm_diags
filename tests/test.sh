# If you get `ModuleNotFoundError: No module named 'acme_diags'`,
# then uncomment the following line and replace <conda-env> with
# the name of your conda environment.
#conda activate <conda-env>
# Use "&&" so the test will fail if either python call fails. 
cd tests/system && python -m unittest && cd .. && python -m unittest && cd ..
