# Precipitation PDF Testing Scripts

This directory contains test scripts for comparing Terai's method vs e3sm_diags implementation.

## Files

### verify_terai.py
- Verifies Terai's pre-calculated PDF files
- Checks model and GPCP PDF statistics
- Quick verification that Terai's PDFs are correct

### test_model_e3sm_diags.py
- Tests e3sm_diags model PDF calculation
- Compares with Terai's model PDF
- Reports differences (should be < 0.1%)

### test_obs_e3sm_diags.py
- Tests e3sm_diags GPCP PDF calculation
- Compares with Terai's GPCP PDF
- **WARNING**: Checks if data is daily vs monthly

## Usage

```bash
cd ~/e3sm_diags/auxiliary_tools/precip_pdf/testing

# Verify Terai's PDFs
python verify_terai.py

# Test model calculation
python test_model_e3sm_diags.py

# Test obs calculation (WARNING: currently using monthly data!)
python test_obs_e3sm_diags.py
```


## Expected Results

- **Model comparison**: < 0.1% difference (numerical precision)
- **Obs comparison**: Should match if using daily GPCP data
