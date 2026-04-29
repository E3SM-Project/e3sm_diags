from __future__ import annotations

import json

from tests.integration.generate_unified_compat_env import (
    build_env_text,
    build_metadata,
    extract_env_dependencies,
    extract_python_version,
    merge_dependency_spec,
    normalize_dependency_spec,
    parse_version_key,
    select_dependency_specs,
    select_latest_nompi_package,
)

CI_ENV_TEXT = """
name: e3sm_diags_ci
channels:
  - conda-forge
dependencies:
  - python >=3.11,<3.14
  - pip
  - setuptools
  - cartopy >=0.17.0,<0.25.0
  - cartopy_offlinedata
  - cf-units
  - dask !=2024.12.0,!=2024.12.1
  - esmf >=8.8.0=nompi_*
  - matplotlib-base >=3.8.2
  - numpy >=2.0.0,<3.0.0
  - pytest
  - pytest-cov
  - uxarray >=2023.3.0
  - xarray >=2024.3.0
  - xcdat >=0.11.1,<1.0.0
  - xesmf >=0.8.7
"""

REPODATA = {
    "info": {"subdir": "linux-64"},
    "packages.conda": {
        "e3sm-unified-1.12.0-nompi_py312_h67c1ab3_0.conda": {
            "name": "e3sm-unified",
            "version": "1.12.0",
            "build": "nompi_py312_h67c1ab3_0",
            "build_number": 0,
            "subdir": "linux-64",
            "timestamp": 1764025560000,
            "depends": [
                "python >=3.12,<3.13.0a0",
                "dask ==2025.9.1",
                "esmf ==8.9.0 nompi_*",
                "esmpy ==8.9.0",
                "matplotlib-base ==3.10.6",
                "numpy >=2.0.0",
                "uxarray >=2024.12.0",
                "xarray ==2025.9.0",
                "xcdat ==0.10.1",
                "xesmf ==0.8.8",
            ],
        },
        "e3sm-unified-1.12.0-nompi_py313_ha8df868_0.conda": {
            "name": "e3sm-unified",
            "version": "1.12.0",
            "build": "nompi_py313_ha8df868_0",
            "build_number": 0,
            "subdir": "linux-64",
            "timestamp": 1764025570000,
            "depends": [
                "python >=3.13,<3.14.0a0",
                "dask ==2025.9.1",
                "esmf ==8.9.0 nompi_*",
                "esmpy ==8.9.0",
                "matplotlib-base ==3.10.6",
                "numpy >=2.0.0",
                "uxarray >=2024.12.0",
                "xarray ==2025.9.0",
                "xcdat ==0.10.1",
                "xesmf ==0.8.8",
            ],
        },
        "e3sm-unified-1.12.0-mpi_openmpi_py313_had67311_0.conda": {
            "name": "e3sm-unified",
            "version": "1.12.0",
            "build": "mpi_openmpi_py313_had67311_0",
            "build_number": 0,
            "subdir": "linux-64",
            "timestamp": 1764025580000,
            "depends": ["python >=3.13,<3.14.0a0"],
        },
        "e3sm-unified-1.11.0rc10-nompi_py310_h1fc0728_0.conda": {
            "name": "e3sm-unified",
            "version": "1.11.0rc10",
            "build": "nompi_py310_h1fc0728_0",
            "build_number": 0,
            "subdir": "linux-64",
            "timestamp": 1739715688759,
            "depends": ["python >=3.10,<3.11.0a0"],
        },
    },
}


class TestGenerateUnifiedCompatEnv:
    def test_extracts_base_env_dependencies(self):
        result = extract_env_dependencies(CI_ENV_TEXT)

        assert "python >=3.11,<3.14" in result
        assert "cf-units" in result
        assert "esmf >=8.8.0=nompi_*" in result
        assert "xesmf >=0.8.7" in result

    def test_normalizes_dependency_spec_to_solver_safe_form(self):
        assert normalize_dependency_spec("python") == "python"
        assert normalize_dependency_spec("numpy >=2.0.0") == "numpy >=2.0.0"
        assert normalize_dependency_spec("esmf ==8.9.0 nompi_*") == "esmf ==8.9.0"
        assert normalize_dependency_spec("xarray ==2025.9.0  # comment") == (
            "xarray ==2025.9.0"
        )

    def test_parses_version_key(self):
        assert parse_version_key("1.12.0") > parse_version_key("1.11.0rc10")
        assert parse_version_key("1.12.0") > parse_version_key("1.12.0rc4")
        assert parse_version_key("1.12.0rc10") > parse_version_key("1.12.0rc4")

    def test_selects_latest_linux_64_nompi_package(self):
        result = select_latest_nompi_package(REPODATA)

        assert result["version"] == "1.12.0"
        assert result["build"] == "nompi_py313_ha8df868_0"
        assert result["filename"] == "e3sm-unified-1.12.0-nompi_py313_ha8df868_0.conda"

    def test_extracts_python_version_from_selected_package(self):
        package = select_latest_nompi_package(REPODATA)

        assert extract_python_version(package) == "3.13"

    def test_merges_base_build_suffix_into_override_spec(self):
        assert (
            merge_dependency_spec("esmf >=8.8.0=nompi_*", "esmf ==8.9.0")
            == "esmf ==8.9.0=nompi_*"
        )
        assert (
            merge_dependency_spec("xcdat >=0.11.1,<1.0.0", "xcdat ==0.10.1")
            == "xcdat ==0.10.1"
        )

    def test_selects_target_dependency_subset(self):
        package = select_latest_nompi_package(REPODATA)

        result = select_dependency_specs(
            package["depends"],
            (
                "python",
                "uxarray",
                "xcdat",
                "dask",
                "esmf",
                "esmpy",
                "matplotlib-base",
                "numpy",
                "xarray",
                "xesmf",
                "pandas",
                "xgcm",
            ),
        )

        assert result == [
            "dask ==2025.9.1",
            "esmf ==8.9.0",
            "esmpy ==8.9.0",
            "matplotlib-base ==3.10.6",
            "numpy >=2.0.0",
            "python >=3.13,<3.14.0a0",
            "uxarray >=2024.12.0",
            "xarray ==2025.9.0",
            "xcdat ==0.10.1",
            "xesmf ==0.8.8",
        ]

    def test_builds_env_text_with_ci_support_packages(self):
        base_dependency_specs = extract_env_dependencies(CI_ENV_TEXT)
        env_text = build_env_text(
            base_dependency_specs=base_dependency_specs,
            dependency_specs=[
                "python >=3.13,<3.14.0a0",
                "esmf ==8.9.0",
                "matplotlib-base ==3.10.6",
                "esmpy ==8.9.0",
                "xarray ==2025.9.0",
            ],
            python_version="3.13",
            env_name="compat",
        )

        assert "name: compat" in env_text
        assert "  - python =3.13" in env_text
        assert "  - pip" in env_text
        assert "  - setuptools" in env_text
        assert "  - pytest" in env_text
        assert "  - pytest-cov" in env_text
        assert "  - cf-units" in env_text
        assert "  - cartopy_offlinedata" in env_text
        assert "  - esmf ==8.9.0=nompi_*" in env_text
        assert "  - esmpy ==8.9.0" in env_text
        assert "  - matplotlib-base ==3.10.6" in env_text
        assert "  - xarray ==2025.9.0" in env_text

    def test_builds_metadata_with_env_hash(self):
        package = select_latest_nompi_package(REPODATA)
        metadata = build_metadata(
            repodata_url="https://conda.anaconda.org/conda-forge/linux-64/repodata.json.bz2",
            package=package,
            dependency_specs=["xarray ==2025.9.0"],
            env_text="name: compat\n",
            python_version="3.13",
        )

        assert metadata["version"] == "1.12.0"
        assert metadata["package_build"] == "nompi_py313_ha8df868_0"
        assert metadata["package_subdir"] == "linux-64"
        assert metadata["python_version"] == "3.13"
        assert metadata["dependency_specs"] == ["xarray ==2025.9.0"]
        assert isinstance(metadata["env_hash"], str)
        assert len(metadata["env_hash"]) == 64

        serialized = json.dumps(metadata, sort_keys=True)
        assert "repodata.json.bz2" in serialized
