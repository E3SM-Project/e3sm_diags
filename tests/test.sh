# Run this script from the top level directory of the repository.

set -e # Fail if any line fails

printf "1. Run unit tests\n"
printf "==============================================\n"
pytest tests/e3sm_diags/

printf "\n2. Run targeted image-regression tests\n"
printf "==============================================\n"
pytest tests/integration/test_plot_image_regressions.py -m image_regression

printf "\n3. Copy over broad integration test data\n"
printf "==============================================\n"
python -m tests.integration.download_data --data-only

printf "\n4. Run broad integration tests\n"
printf "==============================================\n"
CHECK_IMAGES=False pytest tests/integration -m 'not image_regression'
