.PHONY: clean clean-test clean-pyc clean-build docs help test test-unit test-integration test-image-regression test-complete test-complete-validate test-complete-compare promote-complete complete-run-ops-init complete-run-ops-token-create complete-run-ops-env-create complete-run-ops-env-update complete-run-ops-env-show complete-run-scron-config complete-run-scron-validate complete-run-scron-install complete-run-scron-show complete-run-scron-remove
.DEFAULT_GOAL := help

define BROWSER_PYSCRIPT
import os, webbrowser, sys

from urllib.request import pathname2url

webbrowser.open("file://" + pathname2url(os.path.abspath(sys.argv[1])))
endef
export BROWSER_PYSCRIPT

define PRINT_HELP_PYSCRIPT
import re, sys

for line in sys.stdin:
	match = re.match(r'^([a-zA-Z_-]+):.*?## (.*)$$', line)
	if match:
		target, help = match.groups()
		print("%-20s %s" % (target, help))
endef
export PRINT_HELP_PYSCRIPT

BROWSER := python -c "$$BROWSER_PYSCRIPT"

# To run these commands: make <COMMAND>
# ==================================================

help:
	@python -c "$$PRINT_HELP_PYSCRIPT" < $(MAKEFILE_LIST)

# Clean local repository
# ----------------------
clean: clean-build clean-pyc clean-test ## remove all build, test, coverage and Python artifacts

clean-build: ## remove build artifacts
	rm -fr build/
	rm -fr conda-build/
	rm -fr dist/
	rm -fr .eggs/
	find . -name '*.egg-info' -exec rm -fr {} +
	find . -name '*.egg' -exec rm -f {} +

clean-pyc: ## remove Python file artifacts
	find . -name '*.pyc' -exec rm -f {} +
	find . -name '*.pyo' -exec rm -f {} +
	find . -name '*~' -exec rm -f {} +
	find . -name '__pycache__' -exec rm -fr {} +

clean-test: ## remove test and coverage artifacts
	rm -fr tests_coverage_reports/
	rm -f .coverage
	rm -fr htmlcov/
	rm -f coverage.xml
	rm -fr .pytest_cache
	rm -rf .mypy_cache

clean-test-integration: ## remove integration test and artifacts
	rm -rf tests/__pycache__
	rm -rf tests/integration/__pycache__
	rm -rf tests/integration/all_sets_results_test
	rm -rf tests/integration/image_check_failures
	rm -rf tests/integration/integration_test_data
	rm -rf tests/integration/integration_test_images

clean-test-int-res: ## remove integration test results and image check failures
	rm -rf tests/integration/all_sets_results_test
	rm -rf tests/integration/image_check_failures

clean-test-int-data:  # remove integration test data and images (expected) -- useful when they are updated
	rm -rf tests/integration/integration_test_data
	rm -rf tests/integration/integration_test_images

# Quality Assurance
# ----------------------
pre-commit:  # run pre-commit quality assurance checks
	pre-commit run --all-files

lint: ## check style ruff
	ruff check --select I --fix
	ruff check --fix

format: ## format code using ruff
	ruff format

test: ## run tests quickly with the default Python and produces code coverage report
	pytest
	$(BROWSER) tests_coverage_reports/htmlcov/index.html

test-unit: ## run the unit test suite
	pytest tests/e3sm_diags

test-integration: ## download data and run broad integration tests
	python -m tests.integration.download_data --data-only
	CHECK_IMAGES=False pytest tests/integration -m 'not image_regression'

test-image-regression: ## run targeted PNG baseline image-regression tests
	pytest tests/integration/test_plot_image_regressions.py -m image_regression

test-complete: ## run the HPC complete diagnostics workflow
	python -m tests.complete_run.run

test-complete-validate: ## run and compare a complete-run candidate with the accepted baseline
	python -m tests.complete_run.validate

test-complete-compare: ## compare complete-run NetCDF and PNG outputs to the accepted baseline; usage: make test-complete-compare RUN_DIR=/path/to/results [BASELINE_DIR=/path/to/baseline]
	@test -n "$(RUN_DIR)" || { echo "Please specify RUN_DIR=/path/to/results" >&2; exit 2; }
	python -m tests.complete_run.compare --dev-dir "$(RUN_DIR)" $(if $(BASELINE_DIR),--baseline-dir "$(BASELINE_DIR)") --write-diff-pngs --write-diff-html

promote-complete: ## promote reviewed results; usage: make promote-complete RUN_DIR=/path/to/results
	@test -n "$(RUN_DIR)" || { echo "Please specify RUN_DIR=/path/to/results" >&2; exit 2; }
	python -m tests.complete_run.baseline promote --run-dir "$(RUN_DIR)" --channel main

complete-run-ops-init: ## create an operations layout; usage: make complete-run-ops-init OPERATIONS_DIR=/absolute/path [BRANCH=main]
	@test -n "$(OPERATIONS_DIR)" || { echo "Please specify OPERATIONS_DIR=/absolute/path" >&2; exit 2; }
	python -m tests.complete_run.scrontab initialize-operations --operations-dir "$(OPERATIONS_DIR)" --repository-url "$(or $(REPOSITORY_URL),https://github.com/E3SM-Project/e3sm_diags.git)" --branch "$(or $(BRANCH),main)"

complete-run-ops-token-create: ## securely create the controller reporting token; usage: make complete-run-ops-token-create [TOKEN_FILE=$$HOME/.config/e3sm_diags/e3sm_diags-token]
	@TOKEN_FILE="$(or $(TOKEN_FILE),$(HOME)/.config/e3sm_diags/e3sm_diags-token)" bash -c 'set -e; token_file="$$TOKEN_FILE"; install -d -m 700 "$$(dirname "$$token_file")"; read -r -s -p "Paste the E3SM Diags token: " token; printf "\n"; (umask 077; printf "%s\n" "$$token" > "$$token_file"); chmod 600 "$$token_file"; unset token'

complete-run-ops-env-create: ## create the persistent operations environment; usage: make complete-run-ops-env-create CONFIG=/absolute/path/controller.env
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	python -m tests.complete_run.scrontab create-controller-env --config "$(CONFIG)"

complete-run-ops-env-update: ## update the persistent operations environment; usage: make complete-run-ops-env-update CONFIG=/absolute/path/controller.env CONFIRM=YES
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	@test "$(CONFIRM)" = "YES" || { echo "Refusing update; specify CONFIRM=YES" >&2; exit 2; }
	python -m tests.complete_run.scrontab update-controller-env --config "$(CONFIG)" --confirm

complete-run-ops-env-show: ## show persistent operations environment metadata; usage: make complete-run-ops-env-show CONFIG=/absolute/path/controller.env
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	python -m tests.complete_run.scrontab show-controller-env --config "$(CONFIG)"

complete-run-scron-config: ## create an external controller config; usage: make complete-run-scron-config CONFIG=/absolute/path/controller.env
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	python -m tests.complete_run.scrontab create-config --config "$(CONFIG)"

complete-run-scron-validate: ## validate a scheduler config; usage: make complete-run-scron-validate CONFIG=/absolute/path/controller.env
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	python -m tests.complete_run.scrontab validate --config "$(CONFIG)"

complete-run-scron-install: ## install the NERSC scrontab; usage: make complete-run-scron-install CONFIG=/absolute/path/controller.env
	@test -n "$(CONFIG)" || { echo "Please specify CONFIG=/absolute/path/controller.env" >&2; exit 2; }
	python -m tests.complete_run.scrontab install --config "$(CONFIG)"

complete-run-scron-show: ## show the installed NERSC complete-run scrontab
	scrontab -l

complete-run-scron-remove: ## remove the NERSC complete-run scrontab; usage: make complete-run-scron-remove CONFIRM=YES
	@test "$(CONFIRM)" = "YES" || { echo "Refusing removal; specify CONFIRM=YES" >&2; exit 2; }
	python -m tests.complete_run.scrontab remove

# Documentation
# ----------------------
docs: ## generate Sphinx HTML documentation, including API docs
	rm -rf docs/generated
	cd docs && make html
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.htm

docs-versioned: ## generate verisoned Sphinx HTML documentation, including API docs
	rm -rf docs/generated
	cd docs && sphinx-multiversion source _build/html
	$(MAKE) -C docs clean
	$(MAKE) -C docs html
	$(BROWSER) docs/_build/html/index.html

# Build
# ----------------------
env: ## create a conda environment; usage: make env NAME=your_env_name
ifndef NAME
	$(error Please specify the environment name with NAME=your_env_name)
endif
	conda env create -f conda-env/dev.yml -n $(NAME) -y

install: clean ## install the package to the active Python's site-packages
	python -m pip install .
