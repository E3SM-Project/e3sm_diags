# Run this script from the top level directory of the repository.

set -e # Fail if any line fails

# 1. Run unit tests
python -m unittest tests/e3sm_diags/*/test_*.py

# 2. Copy over test data and images from the web. [Takes about four minutes]
# `-m` executes the the module as a stand-alone within a package, allowing for absolute imports
# of other modules within the tests directory.
# https://stackoverflow.com/a/3617928
python -m tests.integration.download_data

# To use `scp`, run the following two lines instead. Be sure to change <username>. [Takes about two minutes]
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/integration_test_data integration_test_data
#scp -r <username>@blues.lcrc.anl.gov:/lcrc/group/e3sm/public_html/e3sm_diags_test_data/integration/expected/integration_test_images integration_test_images

# 3. Run integration tests
python -m unittest discover -s "tests/integration"
