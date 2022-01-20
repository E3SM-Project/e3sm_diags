# Run E3SM Diags on CMIP6 output (on ACME1 which host the complete CMIP archive)

1. Generate xml files

python generate_xmls.py

Run once for ‘amip’, once for ‘historical’. This will generate xml files for each model, one per variable. Only the first realization (r1i1p1f1) is selected. This could be generalized. 
2. Run E3SM Diags

python run_e3sm_diags.py

Run once for ‘amip’, and a second time for ‘historical’. Will submit script to run E3SM Diags on all model output found in step 1 above.

3. Generate webpage

python generate_page.py
