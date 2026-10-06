.PHONY: clean clean-test clean-pyc clean-build compute-node docs help test test-unit test-integration test-image-regression refresh-image-regression test-complete test-complete-validate test-complete-compare promote-complete ops ops-logs ops-update ops-run ops-report ops-help ops-init ops-env ops-enable ops-disable ops-token-create ops-shortcut
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
MACHINE ?= perlmutter

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

# Compute Resources
# -----------------
compute-node: ## request an interactive node; usage: make compute-node MACHINE={anvil,chrysalis,compy,perlmutter} TIME=HH:MM:SS
	@test -n "$(TIME)" || { echo "Please specify TIME=HH:MM:SS" >&2; exit 2; }
	@case "$(MACHINE)" in \
		anvil|chrysalis) srun --pty --nodes=1 --time=$(TIME) /bin/bash ;; \
		compy) salloc --nodes=1 --account=e3sm --time=$(TIME) ;; \
		perlmutter) salloc --nodes 1 --qos interactive --time $(TIME) --constraint cpu --account=e3sm ;; \
		*) echo "Unsupported MACHINE: $(MACHINE). Choose anvil, chrysalis, compy, or perlmutter." >&2; exit 2 ;; \
	esac

# Testing
# -------
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

refresh-image-regression: ## refresh targeted PNG image-regression baselines and metadata
	python -m tests.integration.refresh_plot_image_baselines

test-complete: ## run the HPC complete diagnostics workflow
	python -m tests.complete_run.run

test-complete-validate: ## run and compare a complete-run candidate with the accepted baseline
	python -m tests.complete_run.validate

test-complete-compare: ## compare complete-run NetCDF and PNG outputs to the accepted baseline; usage: make test-complete-compare RUN_DIR=/path/to/results [BASELINE_DIR=/path/to/baseline]
	@test -n "$(RUN_DIR)" || { echo "Please specify RUN_DIR=/path/to/results" >&2; exit 2; }
	python -m tests.complete_run.compare --dev-dir "$(RUN_DIR)" $(if $(BASELINE_DIR),--baseline-dir "$(BASELINE_DIR)") --write-diff-pngs --write-diff-html

promote-complete: ## promote reviewed results; usage: make promote-complete RUN_DIR=/path/to/results [ALLOW_NON_MAIN=1 (maintainer override)]
	@test -n "$(RUN_DIR)" || { echo "Please specify RUN_DIR=/path/to/results" >&2; exit 2; }
	python -m tests.complete_run.baseline promote --run-dir "$(RUN_DIR)" --channel main $(if $(filter 1,$(ALLOW_NON_MAIN)),--allow-non-main)

# Operations
# ----------
# Pass explicit Make inputs through namespaced variables, not shell code.
# Ambient CONFIRM/CONFIG/LINES must not authorize actions or select deployments.
override OPS_INPUT_CONFIG := $(if $(filter command line,$(origin CONFIG)),$(CONFIG))
override OPS_INPUT_CONFIRM := $(if $(filter command line,$(origin CONFIRM)),$(CONFIRM))
override OPS_INPUT_ACTION := $(if $(filter command line,$(origin ACTION)),$(ACTION))
override OPS_INPUT_JOB := $(if $(filter command line,$(origin JOB)),$(JOB))
override OPS_INPUT_LINES := $(if $(filter command line,$(origin LINES)),$(LINES))
override OPS_INPUT_OPERATIONS_DIR := $(if $(filter command line,$(origin OPERATIONS_DIR)),$(OPERATIONS_DIR))
override OPS_INPUT_REPOSITORY_URL := $(if $(filter command line,$(origin REPOSITORY_URL)),$(REPOSITORY_URL))
override OPS_INPUT_BRANCH := $(if $(filter command line,$(origin BRANCH)),$(BRANCH))
override OPS_INPUT_TOKEN_FILE := $(if $(filter command line,$(origin TOKEN_FILE)),$(TOKEN_FILE))
export OPS_INPUT_CONFIG OPS_INPUT_CONFIRM OPS_INPUT_ACTION OPS_INPUT_JOB OPS_INPUT_LINES
export OPS_INPUT_OPERATIONS_DIR OPS_INPUT_REPOSITORY_URL OPS_INPUT_BRANCH OPS_INPUT_TOKEN_FILE

ops: ## read-only operations dashboard
	@python -m tests.complete_run.ops status

ops-logs: ## show recent controller/reporter logs; optional JOB=<id> LINES=100
ops-update: ## fast-forward the operations checkout and validate; no Conda update
ops-run: ## submit a catch-up run; requires CONFIRM=YES
ops-report: ## process completed runs and retry publication; requires CONFIRM=YES
ops-help: ## show everyday operations and administration commands
ops-init: ## create operations layout without scheduling; machine defaults or OPERATIONS_DIR=/path
ops-env: ## explicitly create/update controller environment; ACTION=create or ACTION=update CONFIRM=YES
ops-enable: ## validate deployment and install managed schedule; requires CONFIRM=YES
ops-disable: ## remove only managed schedule, retaining jobs/results; requires CONFIRM=YES
ops-token-create: ## interactively create a private reporting token
ops-shortcut: ## print an optional safely quoted Bash function

ops-logs ops-update ops-run ops-report ops-help ops-init ops-env ops-enable ops-disable ops-token-create ops-shortcut:
	@python -m tests.complete_run.ops $(patsubst ops-%,%,$@)

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
